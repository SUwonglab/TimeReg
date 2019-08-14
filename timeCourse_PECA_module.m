function [TFCluster1,TGCluster1,DriverTF,match]=timeCourse_PECA_module(TRS_norm,Exp,Exp2,Symbol,Sample,TFTG_corr_public,TFTG_corr_RA,lambda,TFName,TFExp,K,Outdir)
R=(1-lambda)*TFTG_corr_public+lambda*TFTG_corr_RA;
m=size(Exp,2);
%%module
for ii=1:m
[TFCluster1{1,ii},TGCluster1{1,ii},W{1,ii},H{1,ii},p1{1,ii},p2{1,ii},FC1{1,ii},FC2{1,ii},ZZ{1,ii}]=TFTG_co_module(TRS_norm{1,ii},R,TFName,Symbol,K(ii),Exp2(:,ii));
K1(ii)=size(TGCluster1{1,ii},2);
if K1(ii)>1
print(1,'-dpng',[Outdir,'/',Sample{ii,1},'_TF_TG_heatmap.png'])
print(2,'-dpng',[Outdir,'/',Sample{ii,1},'_TF_Specific.png'])
close all 
for j=1:K1(ii)
    filename=[Outdir,'/',Sample{ii},'_module',int2str(j),'_Target.txt'];
    fid=fopen(filename,'wt');
    for iter=1:size(TGCluster1{1,ii}{1,j},1)
          fprintf(fid, '%s\t',TGCluster1{1,ii}{1,j}{iter,1});
          fprintf(fid, '%g\t',TGCluster1{1,ii}{1,j}{iter,2});
          fprintf(fid, '%g\t',TGCluster1{1,ii}{1,j}{iter,3});
          fprintf(fid, '%g\t',TGCluster1{1,ii}{1,j}{iter,4});
          fprintf(fid, '%g\n',TGCluster1{1,ii}{1,j}{iter,5});
    end
    fclose(fid);
    filename=[Outdir,'/',Sample{ii},'_module',int2str(j),'_TF.txt'];
    fid=fopen(filename,'wt');
    for iter=1:size(TFCluster1{1,ii}{1,j},1)
          fprintf(fid, '%s\t',TFCluster1{1,ii}{1,j}{iter,1});
          fprintf(fid, '%g\t',TFCluster1{1,ii}{1,j}{iter,2});
          fprintf(fid, '%g\t',TFCluster1{1,ii}{1,j}{iter,3});
          fprintf(fid, '%g\t',TFCluster1{1,ii}{1,j}{iter,4});
          fprintf(fid, '%g\n',TFCluster1{1,ii}{1,j}{iter,5});
    end
    fclose(fid);
end
end
end
%%driver
for ii=2:m
gene1=Symbol(Exp(:,ii)./(Exp(:,ii-1)+0.5)>2);
for j=1:K(ii)
    gene=intersect(gene1,TGCluster1{1,ii}{1,j}(:,1));
    if size(gene,1)>50
        d=ismember(TFName,TFCluster1{1,ii}{1,j}(:,1));
        S2=TRS_norm{1,ii}(d==1,:)-log2(1+TFExp(d==1,ii))-log2(1+Exp(:,ii)');
        [score q1 q3 q4 FC]=Driver(gene,S2,R(d==1,:),TFName(d==1,:),Symbol,TFExp(d==1,ii-1),TFExp(d==1,ii));
         sTFName=TFName(d==1);
        FC=log((1+TFExp(d==1,ii))./(1+TFExp(d==1,ii-1)));
        p=q4;
        p(q1>0.01)=1;
        p(q3>0.01)=1;
        score1=FC.*(-log10(p)).*(FC>log(1.5)).*(p<0.01);
        [d f]=sort(score1,'descend');
        DriverTF{1,ii}{1,j}=[sTFName(f(d>0)) num2cell(q1(f(d>0))) num2cell(q3(f(d>0))) num2cell(p(f(d>0)))];
        else
        DriverTF{1,ii}{1,j}=[];
    end
    if size(DriverTF{1,ii}{1,j},1)>0
        filename=[Outdir,'/',Sample{ii},'_module',int2str(j),'_DriverTF.txt'];
        fid=fopen(filename,'wt');
        for iter=1:size(DriverTF{1,ii}{1,j},1)
              fprintf(fid, '%s\t',DriverTF{1,ii}{1,j}{iter,1});
              fprintf(fid, '%g\t',DriverTF{1,ii}{1,j}{iter,2});
              fprintf(fid, '%g\t',DriverTF{1,ii}{1,j}{iter,3});
              fprintf(fid, '%g\n',DriverTF{1,ii}{1,j}{iter,4});
        end
        fclose(fid);
    end
end
end
%%%matching
match=TimeCourseMatching(DriverTF,TFCluster1,TRS_norm,Symbol,TFName,K);
filename=[Outdir,'/','TimeCourse_ancestor-descendant_mapping.txt'];
fid=fopen(filename,'wt');
for iter=1:size(match,1)
              fprintf(fid, '%s\t',[Sample{match(iter,1)-1},'_module',int2str(match(iter,3))]);
              fprintf(fid, '%s\n',[Sample{match(iter,1)},'_module',int2str(match(iter,2))]);
end
fclose(fid);