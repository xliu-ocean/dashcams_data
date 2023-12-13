function [Rxy,mux,s2x,muy,s2y,k,Nk] = xcovar(x,y,kkk) 


 mux = mean(x,'omitnan');
 muy = mean(y,'omitnan');
 s2x = sum((x-mux).^2,'omitnan');
 s2y = sum((y-muy).^2,'omitnan');

 for ik = 1:100
     Rxy(ik) = sum((x(1:end-ik+1)-mux).*(y(ik:end)-muy),'omitnan');
 end

 k=(1:kkk)-1;
 Nk=kkk;

 return