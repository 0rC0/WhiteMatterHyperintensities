% run LST_lpa with all FLAIRs
% requires glob.m
flairs = glob('/media/drive_s/AG/AG-Floeel-Imaging/02-User/AndreaDellOrco/WhiteMatterLesions/SMARTAGE_BIDS_test/flaircoreg/sub-SMART*/anat/sub-SMART*_FLAIR_roWarped.nii.gz')

for i=1:length(flairs)
disp(flairs{i})
ps_LST_lpa(string(flairs{i}),'',1)
end