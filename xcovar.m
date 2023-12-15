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



%%% tested by a example:
% step 1:
% xtt = 1:400;
% y1t = 0.3*cos(xtt/10)+20;
% y2t = 0.5*cos((xtt-8)/10)+20;
% step 2:
% get lag after using function xcovar
% y2tmp = y2t(lag+1:end);
% y1tmp = y1t(1:end-lag);
% xtttmp = xtt(1:end-lag);
% figure, plot(xtttmp,y1tmp,xtttmp,y2tmp)