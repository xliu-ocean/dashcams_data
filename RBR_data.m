%  Jim Lerczak
%  06 September 2019
%
%  Read in and process RBR CTD data
%


clear all
more off

trbr = readtable('065816_20190819_2038_DataTable.xlsx') ;
N = height(trbr) ;
for ii = 1:N
    dn(ii) = datenum(trbr{ii,1}) ;
    c(ii) = trbr{ii,2} ;
    t(ii) = trbr{ii,3} ;
    p(ii) = trbr{ii,4} ;
    turb(ii) = trbr{ii,5} ;
    sp(ii) = trbr{ii,6} ;
    D(ii) = trbr{ii,7} ;
    s(ii) = trbr{ii,8} ;
    sc(ii) = trbr{ii,9} ;
    da(ii) = trbr{ii,10} ;
    ss(ii) = trbr{ii,11} ;
end

save RBR_data.mat dn c t p turb sp D s sc da ss

%  in this data set, there are three casts.  I plot the here:

p0 = nanmean(p(1:800)) ;
nn = [  805  1019
       3949  4983
      10616 10833] ;
  
figure(1)
clf
co = 'kbrgmc' ;
for ip = 1:3 ;
    mm = nn(ip,1):nn(ip,2) ;
    subplot(1,2,1) ;
    plot(t(mm),-(p(mm)-p0),co(ip)) ;
    hold on
    subplot(1,2,2) ;
    plot(s(mm),-(p(mm)-p0),co(ip)) ;
    hold on
end
    
    

