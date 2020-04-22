
%%%%%%%%%%%%
A=dlmread(Expfile,'\t',0,1);
fid = fopen(Expfile);
C = textscan(fid, '%s %*[^\n]');
fclose(fid);
Symbol=C{1,1};
fileID = fopen(sample_trs_files);
C = textscan(fileID,'%s %s');
fclose(fileID);
Sample=C{1,1};
filenames=C{1,2};
if pre_process_required==1
[Exp,TFName,TGName]=pre_processing(A,Symbol,Sample,filenames,TFName_file,TGName_file,Outdir,species,PECA_Module_dir);
gene=TGName;
else
mkdir(Outdir)
TFName=importdata(TFName_file);
TGName=importdata(TGName_file);
gene=Symbol;
Exp=A;
end
%%%%%%%%%%
m=size(Sample,1);
prior=load ([PECA_Module_dir,'/Prior/TFTG_corr_',species,'.mat']);
[d1 f1]=ismember(TFName,prior.TFName);
[d2 f2]=ismember(TGName,prior.List);
TFTG_corr_public=prior.R2(f1(d1==1),f2(d2==1));
TFName=TFName(d1==1);
TGName=TGName(d2==1);
[d f]=ismember(TFName,gene);
TFExp=Exp(f,:);
[d f]=ismember(TGName,gene);
Exp=Exp(f,:);
Exp2=log2(1+Exp)-repmat(prior.Exp_median(f2(d2)),1,m);
TFTG_corr_private=corr(TFExp',Exp');
TFTG_corr_private(isnan(TFTG_corr_private))=0;
if pre_process_required==1
for i=1:m
    TRS_norm{1,i}=dlmread([Outdir,'/',Sample{i,1},'_TRS.txt'],'\t',0,0);
end
else
for i=1:m
    TRS_norm{1,i}=dlmread(filenames{i,1},'\t',1,1);
end
end
%%%%%%%%%%%%%%%%%%%%%
[TFCluster1,TGCluster1,DriverTF,match]=timeCourse_PECA_module(TRS_norm,Exp,Exp2,TGName,Sample,TFTG_corr_public,TFTG_corr_private,lambda,TFName,TFExp,K,Outdir);
