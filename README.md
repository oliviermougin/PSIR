# PSIR conversion program

You can either use the jupyter notebook or the python script to convert PSIR data acquired from Philips scanners into T1w scans that can be used for freesurfer.

The program takes in the folder where the data is, the prefix of the scan (up until _\_modulus_cphase0X_ [if ptoa] or _\_ph_t732_ [if dcm2niix]) and the suffix of the image (supported formats are .img, .nii and .nii.gz). 

Please note that the program expect the data to have been converted using using 'ptoa -f' or dcm2niix. If another conversion software was used, the names of the files may need to differ (let me know and I will amend it).

# T1map using the PSIR and LUT (including the conversion program)

You can either use the jupyter notebook or the python script to compute a T1map from the PSIR data acquired from Philips scanners. The current protocol includes two protocols dubbed either UK7T (ref) or a newly protocol dubbed Pitts (ref2).

The program takes in the folder where the data is, the prefix of the scan (up until _\_modulus_cphase0X_ [if ptoa] or _\_ph_t732_ [if dcm2niix]) and the suffix of the image (supported formats are .img, .nii and .nii.gz). It also needs the B1map (it only support ptoa at the moment, as the format from dicom still need to be debugged). If the LUT file exists (called LUT_XX.pkl), then it will not be regenerated as it takes one hour or more on a standard computer. The output file will have T1map appended at the end.

Please note that the program expect the data to have been converted using using 'ptoa -f' or dcm2niix. If another conversion software was used, the names of the files may need to differ (let me know and I will amend it).
