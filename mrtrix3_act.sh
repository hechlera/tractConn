#!/bin/bash

### author:  ahechler
### date:    200407

### notes
# - check for correct LUT table
# - -eddy options --slm=linear for 800 shells?
# - interaction with xnat?

#######################################################
######################## TO DO ########################
### - more flexibility with flags for pipeline options
###     - add option to keep/erase previous files
###     - option to run parts of pipeline
### - more error checks and helpful error calls?
### - include set -e again
### - add verbose option with output
### - deringing and ROBLEX brain ext?
###     - move BET to beginning of script to allow early inspection?
###     - visual output into terminal possible? File? 
### - SIFT2?

Usage() {
    cat <<EOF

Usage: mrtrix3_act.sh *args

  Performs all necessary steps for anatomically constrained tractography (ACT) modified after Smith et al. (2013; NeuroImg).
  All functions can be looked up in the MRtrix3 docs: https://mrtrix.readthedocs.io/en/latest/
  The pipeline intersects with various software: MRtrix3, FSL, freesurfer, ANTs
  Data has to be in a BIDS compatible format, with the following structure:
		- \input\sub-{SUBJ_ID}\anat...dwi...
		- \input\data\ including parcellation file (i.e. aparc aseg from freesurfer) or atlas template and look up table
		- only ONE unprocessed dwi and T1 file expected
		
  All steps from inital preprocessing to tracogram (.tck) and connectome (.csv) output are included.
		- .nii -> .mif conversion (MRtrix3 format with extended headers)
		- dwi preprocessing: denoising, eddy, topup, bias correction
		- Constrained Spherical Deconvolution (CSD) with response function estimation
		- Parcellation OR atlas registration: Freesurfer recon-all OR MNI atlas -> individual T1 (i.e. Yeo400)
		- ACT: Registration T1->DWI
		- Tractography: Probabilistic FOD method
		- connectome output as matrix in csv

  Necessary flags:
		-s: subject list as txt or csv with one column specifying the subject IDs (i.e. s014\ns015 for sub-s014 folders)
		-i: input directory with BIDS data and /data subfolder for parcellation files
		-o: output directory
EOF
    exit 1
}

##############################################################################
### ADD FLAGS
#include:
#- start at processing step X -> skip to that part 
#-BET fractional intensity value
#-ACT: streamlines + SIFT streamlines
##############################################################################

# echo usage when called without arguments
[ "$#" == 0 ] && Usage

while getopts "s:i:o:h" opt; do

	case ${opt} in
		
		s)
			SUBJ_LIST=${OPTARG};;
		i)
			PROJ_FOLDER=${OPTARG};;
		o)
			OUT_FOLDER=${OPTARG};;
		h)
			Usage;;
		?)
			Usage;;
			
	esac
	
done			

##############################################################################
### PART 1: ESTIMATE WHITE MATTER RESPONSE FUNCTIONS FOR ALL SUBJECTS
##############################################################################

##############################################################################
### STEP 0: SETUP AND INITIATION

# create output folder for CSD response functions 
RFE_FOLDER=${OUT_FOLDER}/RFE_FOLDER

if [ ! -d ${RFE_FOLDER} ]; then
	mkdir ${RFE_FOLDER}
fi

echo '#############################################################################'
echo 'TractConn: Subject list read, checking for completeness of folders and files'

# iterate through all subjects and sessions based on provided file
while IFS=',' read SUBJ_ID; do

	##############################################################################
	### STEP 1: CHECK AND COPY FILES TO OUTPUT FOLDER
	
	# Set a flag that activates on errors
	FLAG=0
	
	echo "Starting with subject ${SUBJ_ID}"
	# Get folderpaths
	SUBJ_PATH=${PROJ_FOLDER}sub-${SUBJ_ID}
	DWI_PATH=${SUBJ_PATH}/dwi
	ANAT_PATH=${SUBJ_PATH}/anat

	if [[ ! -d ${SUBJ_PATH} || ! -d ${DWI_PATH} || ! -d ${ANAT_PATH} ]]; then
		echo "ERROR: BIDS compatible folderstructure not found. Please check layout"
		FLAG=1
	fi
	
	# Check for atlas files and MNI template
	if [ $(find ${PROJ_FOLDER}data/ -name "Schaefer*.nii*" | wc -l) -ne 1 ]; then
		echo "ERROR: Exactly one atlas from the Schaefer library must be provided"
		FLAG=1
	else
		ATLAS=$(find ${PROJ_FOLDER} -name "Schaefer*.nii*")
	fi
	
	MNI2mm=$(find ${FSLDIR}/data/standard -name "MNI152_T1_1mm_brain.nii.gz")
	
	if [[ -z ${ATLAS} || -z ${MNI2mm} ]]; then
		echo 'ERROR: Yeo atlas and/or MNI2mm template were not found.'
		FLAG=1
	fi
	
	if [ $(find ${ANAT_PATH} -name "*T1w*.nii.gz" | wc -l) -ne 1 ]; then
		echo "Either no T1 or multiple found. Please check input folder for correct naming: *T1w*.nii.gz."
		FLAG=1
	fi
	
	# Search for dwi files
	if [ $(find ${DWI_PATH} -name "*.nii.gz" | wc -l) -gt 1 ]; then
		echo "ERROR: Multiple dwi nifti files found. Only one can be provided"
		FLAG=1
	fi
	
	DWI_FILE=$(find ${DWI_PATH} -name "*.nii.gz") 
	BVEC_FILE=$(find ${DWI_PATH} -name "*.bvec") 
	BVAL_FILE=$(find ${DWI_PATH} -name "*.bval")
	JSON_FILE=$(find ${DWI_PATH} -name "*.json")
	
	if [[ -z ${DWI_FILE} || -z ${BVEC_FILE} || -z ${BVAL_FILE} || -z ${JSON_FILE} ]]; then
		echo "ERROR: DWI data not complete. Original nii, bvec, bval ans json must be available"
		FLAG=1
	fi
	
	if [ ${FLAG} -eq 1 ]; then
		echo "Inital data not complete. Please refer to individual error messages before. Exiting pipeline."
		exit 1
	else
		echo "Initial data complete. Continuing with preprocessing."
	fi
	
	# Check for completeness of output folder
	if [ ! -d ${OUT_FOLDER}/sub-${SUBJ_ID} ]; then
		echo "Creating subject specific output folders"
		mkdir ${OUT_FOLDER}/sub-${SUBJ_ID}
		mkdir ${OUT_FOLDER}/sub-${SUBJ_ID}/anat
		mkdir ${OUT_FOLDER}/sub-${SUBJ_ID}/dwi
	else
		if [ ! -d ${OUT_FOLDER}/sub-${SUBJ_ID}/anat ]; then
			echo "Creating output folder for anatomical files"
			mkdir ${OUT_FOLDER}/sub-${SUBJ_ID}/anat
		fi
		
		if [ ! -d ${OUT_FOLDER}/sub-${SUBJ_ID}/dwi ]; then
			echo "Creating output folder for anatomical files"
			mkdir ${OUT_FOLDER}/sub-${SUBJ_ID}/dwi
		fi
	fi
	
	# Copy T1 file to output folder, change path variable

	echo "Copying T1 over to output folder for further processing"
	MPRAGE_FILE=$(find ${ANAT_PATH} -name "*T1w*.nii.gz")
	cp ${MPRAGE_FILE} ${OUT_FOLDER}/sub-${SUBJ_ID}/anat/
	ANAT_PATH=${OUT_FOLDER}/sub-${SUBJ_ID}/anat
		
	# Copy dwi file to output folder, change path variable
	cp ${DWI_FILE} ${OUT_FOLDER}/sub-${SUBJ_ID}/dwi/
	DWI_FILE=$(find ${OUT_FOLDER}/sub-${SUBJ_ID}/dwi/ -name "*.nii.gz")
	
	# Change working directory
	cd {OUT_FOLDER}
	
	# Added -json_import to store all info in mif header
	mrconvert ${DWI_FILE} ${OUT_FOLDER}/sub-${SUBJ_ID}/dwi/sub-${SUBJ_ID}_dwi.mif -fslgrad ${BVEC_FILE} ${BVAL_FILE} -json_import ${JSON_FILE}
	
	DWI_PATH=${OUT_FOLDER}/sub-${SUBJ_ID}/dwi

	##############################################################################
	### STEP 2: PREPROCESSING
	
	if [ $(find ${DWI_PATH} -name "*dwi.mif" | wc -l) -ne 1 ]; then
		echo "DWI file not found in output/dwi folder. Please check."
		exit 1
	else
		DWI_FILE=$(find ${DWI_PATH} -name "*dwi.mif")
	fi

	# define initial output file names
	DWI_DENOISE="${DWI_FILE%.mif}_dn.mif"
	DWI_NOISEMAP="${DWI_FILE%.mif}_noise.mif"
	
	# denoise main dwi file, get output, name next step
	dwidenoise ${DWI_FILE} ${DWI_DENOISE} -noise ${DWI_NOISEMAP}
	DWI_DN=$(find ${DWI_PATH} -name "*dwi_dn.mif")
	DWI_PREPROC="${DWI_DN%.mif}-preproc.mif"

	# check if reverse encoded dwi file is available
	# pseudocode: if -reverse_enc_flag == y do:
		
		# input and output files
		#DWI_INV=$(find ${DWI_PATH} -name "*dwi0.mif")
		#DWI_DENOISE_INV="${DWI_INV%.mif}_dn.mif"
		#DWI_NOISEMAP_INV="${DWI_INV%.mif}_noise.mif"
		
		# denoising
		#dwidenoise ${DWI_INV} ${DWI_DENOISE_INV} -noise ${DWI_NOISEMAP_INV}
		
		## rpe_ call optimized for an additional b0.mif in reverse PhaseDir
		## -align_seepi concatenates the first b0 of the DWI set with the additional images for topup
	
		#DWI_DN_INV=$(find ${DWI_PATH} -name "*dwi0_dn.mif")
		#dwipreproc ${DWI_DN} ${DWI_PREPROC} -rpe_header -se_epi ${DWI_DN_INV} -align_seepi
		
	# pseudocode: else do:
	
	dwipreproc ${DWI_DN} ${DWI_PREPROC} -rpe_header

	# biascorrect with ANTs
	DWI_PP=$(find ${DWI_PATH} -name "*dwi_dn-preproc.mif") 
	DWI_BCOR="${DWI_PP%.mif}-bcor.mif"
	dwibiascorrect ${DWI_PP} ${DWI_BCOR} -ants
	
	if [[ -z ${DWI_BCOR} ]]; then
		echo "ERROR: Pre-processing could not be completed"
		exit 1
		
	else
		echo "Pre-processing completed. Continuing with response function estimation and CSD"
	fi

	##############################################################################
	### STEP 3: ESTIMATE RESPONSE FUNCTIONS
	
	# Estimate response function and move to group folder
	RFE="${DWI_FILE%.mif}_wm-rfe.txt"
	dwi2response tournier ${DWI_BCOR} ${RFE}
	mv ${RFE} ${RFE_FOLDER}/
	
done < $SUBJ_LIST

##############################################################################
### PART 2: CSD, ACT, connectome creation
##############################################################################

while IFS=',' read SUBJECT; do

	##############################################################################
	### STEP 4: CONSTRAINED SPHERICAL DECONVOLUTION
	
	# calculate group-level response function
	OUT_ARFE="avgwm-rfe.txt"
	cd ${RFE_FOLDER}
	average_response ./*.txt ${OUT_ARFE}	
	
	# get and dilate brain mask to ensure that the mask doesn't cut into relevant tissue for streamline termination
	DWI_BCOR=$(find ${DWI_PATH} -name "*-bcor.mif")
	DWI_MASK="${DWI_FILE%.mif}-preproc_mask.mif"
	dwi2mask ${DWI_BCOR} ${DWI_MASK}

	DWI_MASK=$(find ${DWI_PATH} -name "*_mask.mif")
	DWI_MASK_D="${DWI_MASK%.mif}-dil2.mif"
	maskfilter ${DWI_MASK} dilate -npass 2 ${DWI_MASK_D}

	# CSD
	RFE_FILE=$(find ${RFE_FOLDER} -name "*avgwm-rfe.txt") 
	CSD="${DWI_FILE%.mif}_fod.mif"
	dwi2fod csd ${DWI_BCOR} ${RFE_FILE} ${CSD} -mask ${DWI_MASK_D}
	### verbose optional: check fiber orientation density image - mrview fod.mif -odf.load_sh fod.mif
	
	if [[ -z ${CSD} || -z ${DWI_MASK_D} ]]; then
		echo "ERROR: CSD could not be completed. Check for average response function file. "
		exit 1
	else
		echo "CSD completed for current subject ${SUBJ_ID}. Continuing with registration steps."
	fi
	
	
	##############################################################################
	### STEP 5: T1 -> DWI registration
	
	# get files, prepare output names; transform DWI to nii for registration with FSL
	MPRAGE_FILE=$(find ${ANAT_PATH} -name "*T1w*.nii.gz")
	
	if [ -z ${MPRAGE_FILE} ]; then
		echo 'Error: Base T1 image not found in output folder.'
		exit 1
	fi
	
	MPRAGE_BET="${ANAT_PATH}/sub-${SUBJ_ID}_T1w_bet.nii.gz"
	DWI_FILE=$(find ${DWI_PATH} -name "*dwi_dn-preproc-bcor.mif")

	mrconvert ${DWI_FILE} "${DWI_FILE%.mif}.nii.gz"

	DWI_FILE=$(find ${DWI_PATH} -name "*dwi_dn-preproc-bcor.nii.gz")
	DWI_3D="${DWI_FILE%.nii.gz}-3D.nii.gz"

	# Create folder for transformation matrices
	if [ ! -d ${ANAT_PATH}/registration_files ]; then
		echo "Creating output folder for registration matrices"
		mkdir ${ANAT_PATH}/registration_files
	else
		echo "WARNING: Folder for registration files already present. Old files should be deleted."
	fi

	REG_FOLDER=${ANAT_PATH}/registration_files

	# Get 3D version of 4D DWI file
	mrmath ${DWI_FILE} mean -axis 3 ${DWI_3D}

	# BET on T1
	### WARNING: Standard BET parameters might be to restrictive, resulting in cut-off grey matter
	### Test different fractional intensity values (0 permissive - 1 restrictive; 0.5 = default)
	echo 'WARNING: FSL BET is performed with fractional intensity value of 0.2. The default of 0.5 often cuts into the cortex. Check images afterwards.'
	robustfov -i ${MPRAGE_FILE} -r ${MPRAGE_FILE}
	bet ${MPRAGE_FILE} ${MPRAGE_BET} -f 0.2 -R
	MPRAGE_BET=$(find ${ANAT_PATH} -name "*T1w_bet.nii.gz")

	# get transform matrix T1->3D DWI
	flirt -in ${MPRAGE_BET} -ref ${DWI_3D} -omat ${REG_FOLDER}/T12DWI.mat -dof 6

	# use mrtrix transformation to keep native voxel resolution
	transformconvert ${REG_FOLDER}/T12DWI.mat ${MPRAGE_BET} ${DWI_3D} flirt_import ${REG_FOLDER}/T12DWI.mrtrix
	mrtransform ${MPRAGE_BET} -linear ${REG_FOLDER}/T12DWI.mrtrix "${MPRAGE_BET%.nii.gz}_regDWI.nii.gz"
	
	MPRAGE_REG=$(find ${ANAT_PATH} -name "*regDWI.nii.gz")
	if [ -z ${MPRAGE_REG} ]; then
		echo 'ERROR: T1 to DWI registration was not completed.'
		exit 1
	else
		echo 'T1 to DWI registration completed. Continuing with parcellation or atlas registration'
	fi
	
	##############################################################################
	### STEP 6: ATLAS REGISTRATION / T1 PARCELLATION

	# get transformmatrix from MNI->Sub, apply to atlas file (default affine reg (dof 12))
	flirt -in ${MNI2mm} -ref ${MPRAGE_REG} -omat ${REG_FOLDER}/MNI2SUB.mat
	flirt -in ${ATLAS} -ref ${MPRAGE_REG} -out "${ANAT_PATH}/Yeo400_reg.nii.gz" -init ${REG_FOLDER}/MNI2SUB.mat -applyxfm -interp nearestneighbour
	
	##############################################################################
	### STEP 7: SEGMENTATION
	
	MPRAGE_SEG="${MPRAGE_REG%.nii.gz}_seg.nii.gz"
	MPRAGE_GMWMI="${MPRAGE_REG%.nii.gz}_seg-gmwmi.nii.gz"

	# Segmentation
	5ttgen fsl ${MPRAGE_REG} ${MPRAGE_SEG} -premasked

	# gray matter-white matter interface mask
	MPRAGE_SEG=$(find ${ANAT_PATH} -name "*seg.nii.gz") 
	5tt2gmwmi ${MPRAGE_SEG} ${MPRAGE_GMWMI}
	
	if [ -z ${MPRAGE_SEG} || -z ${MPRAGE_GMWMI} ];then
		echo 'ERROR: Segmentation could not be completed'
		exit 1
	else
		echo 'Segmentation completed. Continuing with ACT and connectome creation.'
	fi
	
	##############################################################################
	### STEP 8: ACT AND CONNECTOME CREATION
	
	# get files (original DWI only for naming), name output
	DWI_FILE=$(find ${DWI_PATH} -name "*dwi.mif")
	# fod file already in ${CSD}
	FOD_FILE=$(find ${DWI_PATH} -maxdepth 1 -name "*dwi_fod.mif")
	# this apparently includes a step before that converts all values to integer. Do testing.
	PARC_FILE=$(find ${ANAT_PATH} -name "Yeo400_reg_int.nii.gz")
	# Seg already in ${MPRAGE_SEG}
	SEG_FILE=$(find ${ANAT_PATH} -name "*T1w*seg.nii.gz")
	# GMWMI already in ${MPRAGE_GMWMI}
	GMWMI_FILE=$(find ${ANAT_PATH} -name "*T1w*seg-gmwmi.nii.gz")

	TCK_OUT="${DWI_FILE%dwi.mif}tracts.tck"
	TCKSFT_OUT="${DWI_FILE%dwi.mif}tracts-sift.tck"
	CONN_OUT="${DWI_FILE%dwi.mif}connectome_zeros.csv"

	# anatomically constrained tractography, seeding at GM-WM interface with 20m streamlines
	tckgen ${CSD} ${TCK_OUT} -act ${MPRAGE_SEG} -seed_gmwmi ${MPRAGE_GMWMI} -backtrack -select 20M

	echo "#########################################################"
	echo "tckgen for ${SUBJ_ID} completed"
	echo "#########################################################"

	# filtering down to 1/4 - 5m streamlines
	tcksift ${TCK_OUT} ${CSD} ${TCKSFT_OUT} -act ${MPRAGE_SEG} -term_number 5m

	echo "#########################################################"
	echo "tcksift for ${SUBJ_ID} completed"
	echo "#########################################################"

	# create connectome matrix
	tck2connectome -symmetric -zero_diagonal ${TCKSFT_OUT} ${PARC_FILE} ${CONN_OUT}

	echo "#########################################################"
	echo "connectome for ${SUBJ_ID} completed. Continuing with next subject"
	echo "#########################################################"
	
done < $SUBJ_LIST
	
