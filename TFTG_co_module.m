function [TFCluster1,TGCluster1,W,H,p1,p2,FC1,FC2,ZZ,Q]=TFTG_co_module(BI,R,TFName,Symbol,k,Exp)
opt = statset('MaxIter',200,'Display','final','TolFun',1e-6);
if k<=1
TFCluster1{1,1}=TFName;
TGCluster1{1,1}=Symbol;
W=[];H=[];p1=[];p2=[];FC1=[];FC2=[];ZZ=[];Q=[];
else
beta=1;
BI_1=log2(1+BI);
Z1=zscore(BI_1);
Z2=zscore(BI_1')';
Z=Z1+Z2;
Z(Z<0)=0;
R(R<0)=0;
ZZ=log2(1+Z).*(R.^beta);
[W,H]=nnmf(ZZ,k,'algorithm','mult','replicates',10,'options',opt);
[d S1]=max(W');
[d S2]=max(H);
BI1=[];
for i=1:size(W,2)
BI1=[BI1;ZZ(S1==i,:)];
end
BI2=[];
for i=1:size(W,2)
BI2=[BI2 BI1(:,S2==i)];
end
BI2(BI2>1)=1;
figure
imagesc(BI2)
colormap(hot)
colorbar
for i=1:size(W,2)
TFCluster{1,i}=TFName(S1==i);
TGCluster{1,i}=Symbol(S2==i);
end
%%%%%%%%%%%%Cluster specific TF,TG
for j=1:k
[h p2(j,:)]=ttest2(ZZ(S1==j,:),ZZ(S1~=j,:),'tail','right');
end
p2=p2';
p2(isnan(p2))=1;
p2(p2<10^(-320))=10^(-320);
for j=1:k
FC2(:,j)=(mean(ZZ(S1==j,:)',2)+0.01)./(mean(ZZ(S1~=j,:)',2)+0.01);
end

ZZ_1=ZZ(:,min(p2')'<0.1/length(Symbol));
S2_1=S2(min(p2')'<0.1/length(Symbol));
for j=1:k
[h p1(j,:)]=ttest2(ZZ_1(:,S2_1==j)',ZZ_1(:,S2_1~=j)','tail','right');
end
p1=p1';
p1(isnan(p1))=1;
p1(p1<10^(-320))=10^(-320);
for j=1:k
FC1(:,j)=(mean(ZZ_1(:,S2_1==j),2)+0.01)./(mean(ZZ_1(:,S2_1~=j),2)+0.01);
end

figure
[d TFidx]=ismember(TFName,Symbol);
size(p1)
size(FC1)
size(Exp(TFidx))
DriverScore=-log10(p1).*FC1.*repmat(Exp(TFidx),1,k);
TGScore=-log10(p2).*FC2.*repmat(Exp,1,k);
for i=1:k
    score=[p2(:,i),FC2(:,i),Exp,TGScore(:,i)];
    [d f]=sort(TGScore(:,i),'descend');
    s=Symbol(f);
    score=score(f,:);
    d=ismember(s,TGCluster{1,i});
    TGCluster1{1,i}=[s(d) num2cell(score(d,:))];
    score=[p1(:,i),FC1(:,i),Exp(TFidx),DriverScore(:,i)];   
    [d f]=sort(DriverScore(:,i),'descend');
    s=TFName(f);
    score=score(f,:);
    d=ismember(s,TFCluster{1,i});
    TFCluster1{1,i}=[s(d) num2cell(score(d,:))];
    subplot(ceil(k/4),4,i)
    plot(0,0)
    for j=1:10
        if length(TFCluster1{1,i})>=j
        text(-1+0.05,1.08-0.19*j,TFCluster1{1,i}(j))
        end
        if length(TFCluster1{1,i})>=j+10
        text(-1+0.75,1.08-0.19*j,TFCluster1{1,i}(10+j))
        end
        if length(TFCluster1{1,i})>=j+20
        text(-1+1.45,1.08-0.19*j,TFCluster1{1,i}(20+j))
        end
    end
    title(['TF cluster',int2str(i),'(',int2str(length(TFCluster1{1,i})),')'])
    set(gca,'XTickLabel',{' '})
    set(gca,'YTickLabel',{' '})
end
set(gcf, 'Position', [0, 0, 400*4 ceil(k/4)*300])
end
