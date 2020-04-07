#!bin/bash

### author: ahechler
### date:	200407

#######################################################
######################## TO DO ########################
### - complete all pipeline steps					###
### - folder structure in XNAT? Do the calls work?  ###
### - subject list file call?						###
###		- second call to the file? Does that work?  ###
### - add flags for pipeline options				###
### 	- add option to keep/erase previous files   ###
### - more error checks and helpful error calls		###
### - add verbose option with output				###
### - deringing and ROBLEX brain ext?				###


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

##############################################################################
### PART 1: ESTIMATE WHITE MATTER RESPONSE FUNCTIONS FOR ALL SUBJECTS      ###
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

	cd ${DWI_PATH}

	DWI_FILE=$(find ${DWI_PATH} -name "*.nii.gz") 
	BVEC_FILE=$(find ${DWI_PATH} -name "*.bvec") 
	BVAL_FILE=$(find ${DWI_PATH} -name "*.bval")
	JSON_FILE=$(find ${DWI_PATH} -name "*.json")

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
### PART 2: CSD, ACT, connectome creation							       ###
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
	
	##############################################################################
	### STEP 6: ATLAS REGISTRATION / T1 PARCELLATION
	
	##############################################################################
	### STEP 7: SEGMENTATION
	
	##############################################################################
	### STEP 7: ACT AND CONNECTOME CREATION
	
done < $SUBJ_LIST
	
	
	
	
	