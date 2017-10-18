# pipeline-descriptor
Wrapper around dogDescriptor.m function

[OUTPUTARGS] = DOGDESCTIPTOR(INPUTARGS)  
  dog = gauss1(sig1)-gauss2(sig1) : difference kernel  
  out = input*dog : convolution  
  out > max(out)/rt : only keep signal > ratio threshold  
  out : [x-y-z-Ifilt-Iraw] : spatial location (0 index) and filter & raw intensity at that location  

Inputs:  
  inputimage: input file can be tif or h5  
  siz: gaussian kernel width  
  sig1: scale of first gaussian  
  sig2: scale of secong gaussian  
  ROI: reject anything outside of BoundingBox  
  rt: threshold ratio  
  withpadding: [1: default]: flag to pre/post pad [siz] size image with
          mirroring or not. Useful to get rid of edge artifacts on the
          image

Outputs:  
  outputfile: text file that has x-y-z-I values as row vectors  
  des: [Nx5]: x-y-z-Ifilt-Iraw row vector  

Examples:  
  dogDescriptor('/nobackup2/mouselight/cluster/2016-09-25/classifier_output/2016-10-02/00/00314/00314-prob.0.h5',...
  '/groups/mousebrainmicro/mousebrainmicro/cluster/Stitching/2016-09-25/Descriptors/13844-prob.0.txt',...
  '[11 11 11]','[3.405500 3.405500 3.405500]','[4.049845 4.049845 4.049845]','[5 1019 5 1531 10 250]','4')
