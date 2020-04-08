#!bin/bash

### author:  ahechler
### date:    200407

#######################################################
######################## TO DO ########################
### - add flags for pipeline options
###     - add option to keep/erase previous files
### - more error checks and helpful error calls
### - include set -e again
### - add verbose option with output
### - folder structure in XNAT? Do the calls work?
### - we need atlas + MNI images somewhere on XNAT
### - subject list file call?
###	    - second call to the file? Does that work?
### - deringing and ROBLEX brain ext?


Usage() {
    cat <<EOF

Usage: mrtrix3_act.sh *args

  Performs all necessary steps for anatomically constrained tractography (ACT) modified after Smith et al. (2013; NeuroImg).
  All functions can be looked up in the MRtrix3 docs: https://mrtrix.readthedocs.io/en/latest/
  The pipeline intersects with various software: MRtrix3, FSL, freesurfer, ANTs
  Data has to be in a BIDS compatible format, with the following structure:
		- \${PROJ_FOLDER}\derivatives\sub-{SUBJ_ID}\ses-{SES_ID}\anat...dwi...
		- only ONE unprocessed dwi and T1 file expected
		
  All steps from inital preprocessing to tracogram (.tck) and connectome (.csv) output are included.
		- .nii -> .mif conversion (MRtrix3 format with extended headers)
		- dwi preprocessing: denoising, eddy, topup, bias correction
		- Constrained Spherical Deconvolution (CSD) with response function estimation
		- Parcellation OR atlas registration: Freesurfer recon-all OR MNI atlas -> individual T1 (i.e. Yeo400)
		- ACT: Registration T1->DWI; T1 
		- Tractography: Probabilistic FOD method
		- connectome output as matrix in csv

  Necessary flags:
		- subject list as csv file (format: 2 columns; 1: subject ID, 2: session ID)
EOF
    exit 1
}

# echo usage when called without arguments
[ "$#" -ne 1 ] && Usage

##############################################################################
### PART 1: ESTIMATE WHITE MATTER RESPONSE FUNCTIONS FOR ALL SUBJECTS
##############################################################################

##############################################################################
### STEP 0: SETUP AND INITIATION

# get base folder. next level must be \sourcedata and \derivatives
PROJ_FOLDER=''

# create output folder for CSD response functions 
RFE_FOLDER=${PROJ_FOLDER}/derivatives/RFE_FOLDER

if [ -d ${PROJECT_DIRECTORY}/derivatives/${RFE_FOLDER} ]; then
	mkdir ${RFE_FOLDER}
fi
	
# get csv list with subjects and sessions
SUBJ_LIST=$1

# iterate through all subjects and sessions based on provided file
while IFS=',' read SUBJECT SESSION; do

	##############################################################################
	### STEP 1: MIF CONVERSION

	# Get folderpaths
	SUBJ_PATH=${PROJ_FOLDER}/sub-${SUBJ_ID}/
	SES_PATH=${SUBJ_PATH}/ses-${SES_ID}/
	DWI_PATH=${SES_PATH}/dwi
	ANAT_PATH=${SES_PATH}/anat

	if [[ ! -d SUBJ_PATH || ! -d SES_PATH || ! -d DWI_PATH || ! -d ANAT_PATH ]]; then
		echo "ERROR: BIDS compatible folderstructure not found. Please check layout"
		exit 1
	fi
		
	cd ${DWI_PATH}

	DWI_FILE=$(find ${DWI_PATH} -name "*.nii.gz") 
	BVEC_FILE=$(find ${DWI_PATH} -name "*.bvec") 
	BVAL_FILE=$(find ${DWI_PATH} -name "*.bval")
	JSON_FILE=$(find ${DWI_PATH} -name "*.json")
	
	if [[ -z DWI_FILE || -z BVEC_FILE || -z BVAL_FILE || -z JSON_FILE ]]; then
		echo "ERROR: DWI data not complete. Original nii, bvec, bval ans json must be available"
		exit 1
	fi

	# Added -json_import to store all info in mif header
	mrconvert ${DWI_FILE} sub-${SUBJ_ID}_ses-${SES_ID}_dwi.mif -fslgrad ${BVEC_FILE} ${BVAL_FILE} -json_import ${JSON_FILE}

	##############################################################################
	### STEP 2: PREPROCESSING

	DWI_FILE=$(find ${DWI_PATH} -name "*dwi.mif")

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

SUBJ_LIST=$1
while IFS=',' read SUBJECT SESSION; do

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
	
	##############################################################################
	### STEP 5: T1 -> DWI registration
	
	# get files, prepare output names; transform DWI to nii for registration with FSL
	MPRAGE_FILE=$(find ${ANAT_PATH} -name "*T1w*.nii.gz")
	MPRAGE_BET="${ANAT_PATH}/sub-${SUBJ_ID}_ses-${SES_ID}_T1w_bet.nii.gz"
	DWI_FILE=$(find ${DWI_PATH} -name "*dwi_dn-preproc-bcor.mif")

	mrconvert ${DWI_FILE} "${DWI_FILE%.mif}.nii.gz"

	DWI_FILE=$(find ${DWI_PATH} -name "*dwi_dn-preproc-bcor.nii.gz")
	DWI_3D="${DWI_FILE%.nii.gz}-3D.nii.gz"

	# Create folder for transformation matrices
	if [ ! -d ${ANAT_PATH}/registration_files ]; then
		echo "Creating output folder for registration matrices"
		mkdir ${ANAT_PATH}/registration_files
	fi

	REG_FOLDER=${ANAT_PATH}/registration_files

	# Get 3D version of 4D DWI file
	mrmath ${DWI_FILE} mean -axis 3 ${DWI_3D}

	# BET on T1
	### WARNING: Standard BET parameters might be to restrictive, resulting in cut-off grey matter
	### Test different fractional intensity values (0 permissive - 1 restrictive; 0.5 = default)
	robustfov -i ${MPRAGE_FILE} -r ${MPRAGE_FILE}
	bet ${MPRAGE_FILE} ${MPRAGE_BET} -f 0.2 -R
	MPRAGE_BET=$(find ${ANAT_PATH} -name "*T1w_bet.nii.gz")

	# get transform matrix T1->3D DWI
	flirt -in ${MPRAGE_BET} -ref ${DWI_3D} -omat ${REG_FOLDER}/T12DWI.mat -dof 6

	# use mrtrix transformation to keep native voxel resolution
	transformconvert ${REG_FOLDER}/T12DWI.mat ${MPRAGE_BET} ${DWI_3D} flirt_import ${REG_FOLDER}/T12DWI.mrtrix
	mrtransform ${MPRAGE_BET} -linear ${REG_FOLDER}/T12DWI.mrtrix "${MPRAGE_BET%.nii.gz}_regDWI.nii.gz"
	
	##############################################################################
	### STEP 6: ATLAS REGISTRATION / T1 PARCELLATION
	
	MPRAGE_REG=$(find ${ANAT_PATH} -name "*regDWI.nii.gz")
	ATLAS=$(find ${PROJECT_DIRECTORY}/derivatives/PARCELLATION -name "Schaefer*.nii.gz")
	MNI2mm=$(find ${PROJECT_DIRECTORY}/derivatives/PARCELLATION -name "MNI152_T1_1mm_brain.nii.gz")

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
	
	##############################################################################
	### STEP 8: ACT AND CONNECTOME CREATION
	
	# get files (original DWI only for naming), name output
	DWI_FILE=$(find ${DWI_PATH} -name "*dwi.mif")
	FOD_FILE=$(find ${DWI_PATH} -maxdepth 1 -name "*dwi_fod.mif")
	PARC_FILE=$(find ${ANAT_PATH} -name "Yeo400_reg_int.nii.gz")
	SEG_FILE=$(find ${ANAT_PATH} -name "*T1w*seg.nii.gz")
	GMWMI_FILE=$(find ${ANAT_PATH} -name "*T1w*seg-gmwmi.nii.gz")

	TCK_OUT="${DWI_FILE%dwi.mif}tracts.tck"
	TCKSFT_OUT="${DWI_FILE%dwi.mif}tracts-sift.tck"
	CONN_OUT="${DWI_FILE%dwi.mif}connectome_zeros.csv"

	# anatomically constrained tractography, seeding at GM-WM interface with 20m streamlines
	tckgen ${FOD_FILE} ${TCK_OUT} -act ${SEG_FILE} -seed_gmwmi ${GMWMI_FILE} -backtrack -select 20M

	echo "#########################################################"
	echo "tckgen for ${SUBJ_ID} completed"
	echo "#########################################################"

	# filtering down to 1/4 - 5m streamlines
	tcksift ${TCK_OUT} ${FOD_FILE} ${TCKSFT_OUT} -act ${SEG_FILE} -term_number 5m

	echo "#########################################################"
	echo "tcksift for ${SUBJ_ID} completed"
	echo "#########################################################"

	# create connectome matrix
	tck2connectome -symmetric -zero_diagonal ${TCKSFT_OUT} ${PARC_FILE} ${CONN_OUT}

	echo "#########################################################"
	echo "connectome for ${SUBJ_ID} completed"
	echo "#########################################################"
	
done < $SUBJ_LIST
	
