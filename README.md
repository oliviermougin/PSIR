# PSIR conversion program

You can either use the jupyter notebook or the python script to convert PSIR data acquired from Philips scanners into T1w scans that can be used for freesurfer.

The program takes in the folder where the data is, the prefix of the scan (up until _\_modulus_cphase0X_ [if ptoa] or _\_ph_t732 [if dcm2niix]) and the suffix of the image (supported formats are .img, .nii and .nii.gz). 

Please note that the program expect the data to have been converted using using 'ptoa -f' or dcm2niix. If another conversion software was used, the names of the files may need to differ (let me know and I will amend it).
