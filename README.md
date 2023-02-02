# BSS_ADJUST

repository with code that process BSS+ADJUST EEG artifact rejection on Matlab for EKSO gait trials //

Be sure you have a verson of EEGlab intalled in your Matlab prompt that is newer than the 14.1.1. You can download and install it from the EEGlab download page [https://sccn.ucsd.edu/eeglab/downloadtoolbox.php/download.php](https://sccn.ucsd.edu/eeglab/downloadtoolbox.php/download.php). After that please locate your eeglab path, and install the following plugin dependencies.

- Cleanline 2.00
- clean_rawdata 2.3
- ADJUST 1.1.1
- Fieldtrip-lite
- AAR 13.11

You can go to the EEG plugins page and donwload them directly [https://sccn.ucsd.edu/eeglab/plugin_uploader/plugin_list_all.php](https://sccn.ucsd.edu/eeglab/plugin_uploader/plugin_list_all.php), or use the eeglab command and the option Manage Plugins in your Matlab prompt. 

For executing the artifact rejection using BSS+ADJUST simply run the following command having your EKSO .mat files added on your Matlab path.

```matlab
reading_gtec_EKSO_mode('/path_to_EKSO_Data/S9_ekso_modeA_2.mat','modeA_2')
```
Or use the general code to calculate ERSP and CMC that used ASR and ADJUST for eliminating artifacts as the following command
