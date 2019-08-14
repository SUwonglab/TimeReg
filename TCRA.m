
Sample=importdata(SampleNamefile);
m=size(Sample,1);
TFName=importdata(TFNamefile);
TGName=importdata(TGNamefile);
prior=load ([PECA_Module_dir,'/Prior/TFTG_corr_',spesies,'.mat']);
[d1 f1]=ismember(TFName,prior.TFName);
[d2 f2]=ismember(TGName,prior.List);
TFTG_corr_public=prior.R2(f1,f2);
Exp=dlmread(Expfile,'\t',0,1);
fid = fopen(Expfile);
C = textscan(fid, '%s %*[^\n]');
fclose(fid);
gene=C{1,1};
[d f]=ismember(TFName,gene);
TFExp=Exp(f,:);
[d f]=ismember(TGName,gene);
Exp=Exp(f,:);
Exp2=log2(1+Exp)-repmat(prior.Exp_median(f2),1,m);
TFTG_corr_private=corr(TFExp',Exp');
for i=1:m
    TRS_norm{1,i}=dlmread([TRS_dir,'/',Sample{i,1},'_TRS.txt']);
end
%%%%%%%%%%%%%%%%%%%%%
[TFCluster1,TGCluster1,DriverTF,match]=timeCourse_PECA_module(TRS_norm,Exp,Exp2,TGName,Sample,TFTG_corr_public,TFTG_corr_private,lambda,TFName,TFExp,K,Outdir);