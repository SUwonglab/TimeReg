function [score p1 p3 p4 FC]=Driver(gene,BI2,R,TFName,Symbol,TFExp1,TFExp2,back)
if nargin < 9
    back=Symbol(rand(length(Symbol),1)<length(gene)/length(Symbol));
    %back=setdiff(Symbol,gene);
end
d=ismember(Symbol,gene);
d1=ismember(Symbol,back);
for i=1:length(TFName)
p1(i,1)=ranksum(R(i,d==1),R(i,d1==1),'tail','right');
p1_FC(i,1)=max(mean(R(i,d==1)),0.01)/max(mean(R(i,d1==1)),0.01);
end
p1(isnan(p1))=1;
p1(p1<10^(-320))=10^(-320);
[h crit_p p1]=fdr_bh(p1);

BI2_rank=quantilenorm(BI2);
for i=1:length(TFName)
p3(i,1)=ranksum(BI2_rank(i,d==1),BI2_rank(i,d1==1),'tail','right');
p3_FC(i,1)=mean(BI2_rank(i,d==1))/mean(BI2_rank(i,d1==1));
end
p3(isnan(p3))=1;
p3(p3<10^(-320))=10^(-320);
[h crit_p p3]=fdr_bh(p3);

R1=R;
R1(R<0)=0;
BI3_rank=quantilenorm(BI2).*R1;
for i=1:length(TFName)
p4(i,1)=ranksum(BI3_rank(i,d==1),BI3_rank(i,d1==1),'tail','right');
p4_FC(i,1)=mean(BI3_rank(i,d==1))/mean(BI3_rank(i,d1==1));
end
p4(isnan(p4))=1;
p4(p4<10^(-320))=10^(-320);
[h crit_p p4]=fdr_bh(p4);

FC=log2(TFExp2+1)-log2(TFExp1+1);
%FC(FC<1)=0;
% score=(-log10(p1)-log10(p2)-log10(p3))+FC;
% score=score.*(p1<0.1).*(p2<0.1).*(p3<0.1).*(FC>1.5);
score=(sqrt(-log10(p1).*p1_FC)+sqrt(-log10(p3).*p3_FC))+FC;
score=score.*(p1<0.05).*(p3<0.05).*(FC>log2(1.5));
