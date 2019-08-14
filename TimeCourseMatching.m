function match=TimeCourseMatching(DriverTF,TFCluster1,TRS_norm,Symbol,TFName,K)
match=[];
for ii=2:length(K)
for j=1:K(ii)
if size(DriverTF{1,ii}{1,j},1)>0
s=[DriverTF{1,ii}{1,j}(:,1);TFCluster1{1,ii}{1,j}(:,1)];
else
s=TFCluster1{1,ii}{1,j}(:,1);
end
s=s(1:20,1);
z1=zscore(TRS_norm{1,ii-1}')';
[d f]=ismember(s,Symbol);
z11=z1(:,f);
clear c
clear c1
for i1=1:K(ii-1)
[d1 f1]=ismember(TFCluster1{1,ii-1}{1,i1}(:,1),TFName);
c{1,i1}=z11(f1(1:20),:);
c1(1,i1)=mean(mean(c{1,i1}));
end
[~,match_idx]=max(c1);
match=[match;[ii j match_idx]];
end
end