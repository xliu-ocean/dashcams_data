function [Rxy,mux,s2x,muy,s2y,k,Nk] = xcovar(x,y,kkk) 


 mux = mean(x,'omitnan');
 muy = mean(y,'omitnan');
 s2x = sum((x-mux).^2,'omitnan');
 s2y = sum((y-muy).^2,'omitnan');

 Rxy = sum((x-mux).*(y-muy),'omitnan');

 k=1:kkk;
 Nk=kkk;

 return