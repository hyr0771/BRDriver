function [final_R0] = lioness_method(data,i)
%construct sample specific network

PCC=corrcoef(data');
v0=data;
v0(:,i)=[];
PCC1=corrcoef(v0');
R=size(data,2)*(PCC-PCC1)+PCC1;

final_R = R;
final_R0 = abs(R);
index =isnan(final_R);
final_R(index(1,:),:)=[];  
final_R(:,index(:,1))=[];
A_fR = abs(final_R(:));
threshold = mean( A_fR );
final_R0(isnan(final_R0))=0;
final_R0(final_R0 < threshold)=0;
final_R0(final_R0~=0)=1;

end

