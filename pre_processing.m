function [Exp,TFName,TGName]=pre_processing(A,Symbol,sample,filenames,TFName_file,TGName_file,out_folder,spesies,PECA_Module_dir)
mkdir(out_folder)
specific_list=specific_gene_list(A,Symbol,spesies,PECA_Module_dir);
specific_list=[sample';specific_list];
%xlswrite([out_folder,'/specific_list.xlsx'],specific_list,'specific_list_rank')
filename=[out_folder,'/specific_list.txt'];
fid=fopen(filename,'wt');
for i=1:size(specific_list,1)
for j=1:size(specific_list,2)-1
fprintf(fid, '%s\t',specific_list{i,j});
end
fprintf(fid, '%s\n',specific_list{i,j+1});
end
fclose(fid);
Symbol_selected=unique(specific_list);
TFName=importdata(TFName_file);
TGName=importdata(TGName_file);
[d1 f1]=ismember(TFName,Symbol_selected);
[d2 f2]=ismember(TGName,Symbol_selected);
TFName=TFName(d1==1);
TGName=TGName(d2==1);
for i=1:size(filenames,1)
    TRS=dlmread(filenames{i,1},'\t',1,1);
    TRS_filtered{1,i}=TRS(d1,d2);
end
[d f]=ismember(Symbol,TGName);
Exp=A(d,:);
%%%write
filename=[out_folder,'/TFName.txt'];
 fid=fopen(filename,'wt');
for i=1:size(TFName,1)
 fprintf(fid, '%s\n',TFName{i,1});
 end
 fclose(fid);
filename=[out_folder,'/TGName.txt'];
 fid=fopen(filename,'wt');
for i=1:size(TGName,1)
 fprintf(fid, '%s\n',TGName{i,1});
 end
 fclose(fid);
for i=1:size(filenames,1)
    filename=[out_folder,'/',sample{i,1},'_TRS.txt'];
    dlmwrite(filename,TRS_filtered{1,i},'\t')
    %T = array2table(TRS_filtered{1,i},'RowNames',TFName);
    %writetable(T,filename,'Delimiter','\t','WriteRowNames',true,'WriteVariableNames',true)
end