function specific_list=specific_gene_list(A,Symbol,spesies,PECA_Module_dir)
cut =3000;
load([PECA_Module_dir,'/Prior/TFTG_corr_',spesies,'.mat'])
geneName=intersect(List,Symbol);
[d f]=ismember(geneName,List);
Exp_median=Exp_median(f);
[d f]=ismember(geneName,Symbol);
A=A(f,:);
A_norm=log2(1+quantilenorm(A));
Exp_median_timeCourse=median(A_norm')';
Exp_fold=A_norm-repmat(sqrt(Exp_median.*Exp_median_timeCourse),1,size(A,2));
specific_list=[];
for i=1:size(A,2)
    [d f]=sort(Exp_fold(:,i),'descend');
    specific_list=[specific_list geneName(f(1:cut))];
end