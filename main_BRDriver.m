clc,clear


Data = 'exam';
coding_tumor =strcat(Data,'_coding_tumor.txt'); 
lncRNA_tumor_profile = strcat(Data,'_lncRNA_tumor.txt');
name = coding_tumor(1:5);

[coding_result,lncRNA_result] = BRDriver( coding_tumor ,lncRNA_tumor_profile);

writecell( coding_result,strcat(name,'result_coding.txt'),'Delimiter','tab');
writecell( lncRNA_result,strcat(name,'result_lncRNA.txt'),'Delimiter','tab');

clear coding_normal coding_tumor lncRNA_normal_profile lncRNA_tumor_profile name result 


