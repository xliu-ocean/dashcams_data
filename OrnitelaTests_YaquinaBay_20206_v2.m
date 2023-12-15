%  Jim Lerczak
%  26 June 2020
%
%  Read in and process RBR CTD data
%
%  Revised by Xiaohui Liu
%  Dec 5 2023
%% load data
clear all
more off

dz = 0.25 ;
Nrm = 5 ;

Tmin = 9.5 ;
Tmax = 22 ;
Smax = 34 ;
Smin = 0 ;
Cmin = 0 ;
Cmax = 40 ;

dnmin = datenum(['25-Jun-2020 23:01:18']) ;
dnmax = datenum(['26-Jun-2020 01:46:17']) ;

%  load the RBR data
load Ornitela_Yaquina_2020_06/RBR_data_YaquinaBay_202006.mat
%  load the Ornitela data
load Ornitela_Yaquina_2020_06/Ornitela_casts/Ornitela_YaquinaBay_062020_UAE.mat
%%
% figure
% clf


for ii = 1:4

    
    Drbr = interp1(dn,D,dn10Hz{ii}) ;
    %  get an average data point offset b/w the bird tag and the RBR CTD
    [Rxy,mux,s2x,muy,s2y,k,Nk] = xcovar(Drbr,D10Hz{ii},100) ;
    rhoxy = Rxy./sqrt(s2x.*s2y) ;
    [vl,idx] = max(rhoxy) ;
    lag(ii) = k(idx) ;
    
%     Dtmp{ii}(:) = D10Hz{ii}(lag(ii)+1:end) ;
%     Ttmp{ii}(:) = T10Hz{ii}(lag(ii)+1:end) ;
%     Ctmp{ii}(:) = C10Hz{ii}(lag(ii)+1:end) ;
    Dtmp = D10Hz{ii}(lag(ii)+1:end) ;
    Ttmp = T10Hz{ii}(lag(ii)+1:end) ;
    Ctmp = C10Hz{ii}(lag(ii)+1:end) ;
    dntmp = dn10Hz{ii}(1:end-lag(ii)) ;
%     Drbr2{ii}(:) = interp1(dn,D,dntmp) ;
    Drbr = interp1(dn,D,dntmp) ;
    
    %  Get a mean depth offset between the bird tag and the RBR CTD    
%     dD(ii) = mean(Drbr2{ii}-Dtmp{ii},'omitnan') ;
    dD(ii) = mean(Drbr-Dtmp,'omitnan') ;

%     plot(Drbr,Dtmp+dD(ii),'ko') ;
%     plot(Drbr,Dtmp,'ko') ;
%     hold on

    % interpolate T and C of RBR to Ornitela

     Trbr{ii}(:) = interp1(dn,T,dntmp) ;
     Crbr{ii}(:) = interp1(dn,C,dntmp) ;


end
% plot([0 12],[0 12],'r') ;
% title('Water Depth (Pressure) Comparison')
% ylabel('Ornitela Tag (meters; all tags combined)') ;
% xlabel('RBR CTD (meters)') ;

%%
figure(1)
% print -djpeg90 -r300 OrnitelaTests_YaquinaBay_202006_Pressure.jpg

figure(2)
clf
% plot_size(1,1,10,12) ;

figure(3)
clf
% plot_size(1,1,10,12) ;

figure(4)
clf
% plot_size(1,5,15,5) ;

co = ['kbrgmckbrgmckbrgmckbrgmckbrgmc'] ;
% Now let's analyze temperature for each downcast for each sensor
for ii = 1:4 
    
    figure(4)
    nngps = find((dngps{ii}>=dnmin).*(dngps{ii}<=dnmax)) ;
    plot(dngps{ii}(nngps),Tint{ii}(nngps),[co(ii) '+']) ;
    hold on
    
    figure(2)
    clf
    figure(3) 
    clf

    Ctag = [] ;
    Ttag = [] ;
    Dtag = [] ;
    Crbr = [] ;
    Trbr = [] ;
    DDrbr = [] ;
    
    Cmaxtag = Cmax ;
            
    %  grab the tag data
    dntag = dn10Hz{ii}(1:end-lag(ii)) ;
    dtag = D10Hz{ii}(lag(ii)+1:end) ;
    ttag = T10Hz{ii}(lag(ii)+1:end) ;
    ctag = C10Hz{ii}(lag(ii)+1:end) ;
    
    %  Get a mean depth offset between the bird tag and the RBR CTD    
    dD(ii) = mean(Drbr-Dtmp,'omitnan') ;

    plot(Drbr,Dtmp+dD(ii),'ko') ;
    hold on
            
    for il = 1:5
        for ic = 1:5 ;
            jj = 5*(il - 1) + ic ;
            nn = nn1(jj):nn2(jj) ;
            dndn = dn(nn) ;
            cc = C(nn) ;
            ss = S(nn) ;
            tt = T(nn) ;
            dd = D(nn) ;
            mm = find(dd>=dz) ;
            dndn = dndn(mm(1):end) ;            
            cc = cc(mm(1):end) ;
            ss = ss(mm(1):end) ;
            tt = tt(mm(1):end) ;
            dd = dd(mm(1):end) ;
            mm = find(dd<=(max(dd)-dz)) ;
            dndn = dndn(mm) ;
            cc = cc(mm) ;
            ss = ss(mm) ;
            tt = tt(mm) ;
            dd = dd(mm) ;
            dndn = movmean(dndn,Nrm) ;
            ss = movmean(ss,Nrm) ;
            tt = movmean(tt,Nrm) ;
            dd = movmean(dd,Nrm) ;            
            
            %  plot only the downcasts:

            [vl,idx] = max(dd) ;
            mm = 1:idx ;
            dndn = dndn(mm) ;
            cc = cc(mm) ;
            ss = ss(mm) ;
            tt = tt(mm) ;
            dd = dd(mm) ;
            
            %  get the tag data
            mm = find((dntag>=dndn(1)).*(dntag<=dndn(end))) ;
            dntmp = dntag(mm) ;
            dtmp = dtag(mm) ;
            ttmp = ttag(mm) ;
            ctmp = ctag(mm) ;
            dntmp = movmean(dntmp',Nrm) ;
            dtmp = movmean(dtmp',Nrm) ;
            ttmp = movmean(ttmp',Nrm) ;
            ctmp = movmean(ctmp',Nrm) ;
            
            Dtag = [Dtag dtmp'] ;
            Ttag = [Ttag ttmp'] ;
            Ctag = [Ctag ctmp'] ;
            DDrbr = [DDrbr interp1(dn,D,dntmp)'] ;
            Trbr = [Trbr interp1(dn,T,dntmp)'] ;
            Crbr = [Crbr interp1(dn,C,dntmp)'] ;
            
            if ii == 1
                figure(4)
                plot([1 1]*mean(dndn,'omitnan'),[0 50],'k') ;
            end

            %  Temperature
            figure(2)
            subplot(5,2,(il - 1)*2 + 1) ;
            plot(tt,-dd,co(ic)) ;
            hold on
            title(['Site ' num2str(il) ' (downdast, RBR)']) ;
            ylabel('Water Depth (m)') ;
            subplot(5,2,(il - 1)*2 + 2) ;
            plot(ttmp,-dtmp,co(ic)) ;
            hold on
            title(['Site ' num2str(il) ' (downdast, tag ' num2str(SN(ii)) ')']) ;

            %  Conductivity
            figure(3)
            subplot(5,2,(il - 1)*2 + 1) ;
            plot(cc,-dd,co(ic)) ;
            hold on
            title(['Site ' num2str(il) ' (downdast, RBR)']) ;
            ylabel('Water Depth (m)') ;
            subplot(5,2,(il - 1)*2 + 2) ;
            plot(ctmp,-dtmp,co(ic)) ;
            hold on
            title(['Site ' num2str(il) ' (downdast, tag ' num2str(SN(ii)) ')']) ;
            
            Cmaxtag = max([Cmaxtag max(ctmp)]) ;
            
        end
        
        figure(2)
        subplot(5,2,(il - 1)*2 + 1) ;
        axis([Tmin Tmax -Dmax(il) 0]); 
        subplot(5,2,(il - 1)*2 + 2) ;
        axis([Tmin Tmax -Dmax(il) 0]); 
        
        figure(3)
        subplot(5,2,(il - 1)*2 + 1) ;
        axis([Cmin Cmax -Dmax(il) 0]); 
        subplot(5,2,(il - 1)*2 + 2) ;
        axis([Cmin Cmaxtag -Dmax(il) 0]); 
        
    
    end
    figure(2)
    subplot(5,2,9)
    xlabel('Temperature (^oC)')
    subplot(5,2,10)
    xlabel('Temperature (^oC)')
    fnT = ['OrnitelaTests_YaquinaBay_202006_Temperature_' num2str(SN(ii)) '.jpg'] ;
    eval(['print -djpeg90 -r300 ' fnT]) ;
    figure(3)
    subplot(5,2,9)
    xlabel('Conductivity (mS/cm)')
    subplot(5,2,10)
    xlabel('Conductivity (mS/cm)')
    fnC = ['OrnitelaTests_YaquinaBay_202006_Conductivity_' num2str(SN(ii)) '.jpg'] ;
    eval(['print -djpeg90 -r300 ' fnC]) ;

end

figure(4)
datetick('x') ;
axis([dnmin dnmax 9 35])
title('Tag Interior Temperature vs. Time (Yaquina Bay, 06/25/2020)') ;
ylabel('T (^oC)') 
xlabel('Time') ;
print  -djpeg90 -r300 OrnitelaTests_YaquinaBay_202006_TagInteriorTemp.jpg

figure(5)
clf
% plot_size(1,5,15,5) ;
figure(6)
clf
% plot_size(1,5,15,5) ;
dn1 = datenum('26-Jun-2020 01:29:11') ;
dn2 = datenum('26-Jun-2020 01:38:30') ;
for ii = 1:4 ;
    nn = find((dn10Hz{ii}>=dn1).*(dn10Hz{ii}<=dn2)); 
    dntmp = dn10Hz{ii}(nn) ;
    Ctmp = C10Hz{ii}(nn) ;
    Ttmp = T10Hz{ii}(nn) ;
    figure(5)
    plot(dntmp,Ttmp,[co(ii) '.']) ;
    hold on
    figure(6)
    plot(dntmp,Ctmp,[co(ii) '.']) ;
    hold on
end

figure(5)
datetick('x') ;
axis([dn1 dn2 20.2 22.5])
title('Temperature vs Time (Yaquina Bay, Site 6, 06/25/2020)','fontsize',18) ;
ylabel('T (^oC)','fontsize',18) 
xlabel('Time','fontsize',18) ;
print  -djpeg90 -r300 OrnitelaTests_YaquinaBay_202006_TvsTime.jpg

figure(6)
datetick('x') ;
axis([dn1 dn2 4 16])
title('Conductivity vs Time (Yaquina Bay, Site 6, 06/25/2020)','fontsize',18) ;
ylabel('C (mS/cm)','fontsize',18) 
xlabel('Time','fontsize',18) ;
print  -djpeg90 -r300 OrnitelaTests_YaquinaBay_202006_CvsTime.jpg


% return

%%

% co = ['kbrgmckbrgmckbrgmckbrgmckbrgmc'] ;
% Now let's analyze temperature for each downcast for each sensor
for ii = 1:4
    
%     figure(4)
    nngps = find((dngps{ii}>=dnmin).*(dngps{ii}<=dnmax)) ;
%     plot(dngps{ii}(nngps),Tint{ii}(nngps),[co(ii) '+']) ;
%     hold on
    
%     figure(2)
%     clf
%     figure(3) 
%     clf

    Ctag = [] ;
    Ttag = [] ;
    Dtag = [] ;
    Crbr = [] ;
    Trbr = [] ;
    DDrbr = [] ;
    
    Cmaxtag = Cmax ;
            
    %  grab the tag data
    dntag = dn10Hz{ii}(1:end-lag(ii)) ;
    dtag = D10Hz{ii}(lag(ii)+1:end) ;
    ttag = T10Hz{ii}(lag(ii)+1:end) ;
    ctag = C10Hz{ii}(lag(ii)+1:end) ;
    
    %  Get a mean depth offset between the bird tag and the RBR CTD    
    dD(ii) = mean(Drbr-Dtmp,'omitnan') ;

    plot(Drbr,Dtmp+dD(ii),'ko') ;
    hold on
            
    for il = 1:5
        for ic = 1:5 ;
            jj = 5*(il - 1) + ic ;
            nn = nn1(jj):nn2(jj) ;
            dndn = dn(nn) ;
            cc = C(nn) ;
            ss = S(nn) ;
            tt = T(nn) ;
            dd = D(nn) ;
            mm = find(dd>=dz) ;
            dndn = dndn(mm(1):end) ;            
            cc = cc(mm(1):end) ;
            ss = ss(mm(1):end) ;
            tt = tt(mm(1):end) ;
            dd = dd(mm(1):end) ;
            mm = find(dd<=(max(dd)-dz)) ;
            dndn = dndn(mm) ;
            cc = cc(mm) ;
            ss = ss(mm) ;
            tt = tt(mm) ;
            dd = dd(mm) ;
            dndn = movmean(dndn,Nrm) ;
            ss = movmean(ss,Nrm) ;
            tt = movmean(tt,Nrm) ;
            dd = movmean(dd,Nrm) ;  
            cc = movmean(cc,Nrm) ;           

            % %%%%%% get both the downcasts(***_a) and upper casts(***_b) :

            [vl,idx] = max(dd) ;
            mm = 1:idx ;
            mm2 = idx:length(dd) ;
%             mm = idx:length(dd) ;
%             mm = 1:length(dd);
            dndn1 = dndn(mm) ;
            cc1 = cc(mm) ;
            ss1 = ss(mm) ;
            tt1 = tt(mm) ;
            dd1 = dd(mm) ;

            cc_a{jj,ii} = cc1;
            dndn_a{jj,ii} = dndn1;
            tt_a{jj,ii} = tt1;
            dd_a{jj,ii} = dd1;
            ss_a{jj,ii} = ss1;

            %   %% get upward casts
            dndn2 = dndn(mm2) ;
            cc2 = cc(mm2) ;
            ss2 = ss(mm2) ;
            tt2 = tt(mm2) ;
            dd2 = dd(mm2) ;

            cc_b{jj,ii} = cc2;
            dndn_b{jj,ii} = dndn2;
            tt_b{jj,ii} = tt2;
            dd_b{jj,ii} = dd2;
            ss_b{jj,ii} = ss2;
            
            %  get the tag data
            mm = find((dntag>=dndn1(1)).*(dntag<=dndn1(end))) ;
            dntmp = dntag(mm) ;
            dtmp = dtag(mm) ;
            ttmp = ttag(mm) ;
            ctmp = ctag(mm) ;
            dntmp = movmean(dntmp',Nrm) ;
            dtmp = movmean(dtmp',Nrm) ;
            ttmp = movmean(ttmp',Nrm) ;
            ctmp = movmean(ctmp',Nrm) ;

            dntmp_a{jj,ii} = dntmp;
            dtmp_a{jj,ii} = dtmp;
            ttmp_a{jj,ii} = ttmp;
            ctmp_a{jj,ii} = ctmp;
            
            %  get the tag data (upward casts)
            mm2 = find((dntag>=dndn2(1)).*(dntag<=dndn2(end))) ;
            dntmp = dntag(mm2) ;
            dtmp = dtag(mm2) ;
            ttmp = ttag(mm2) ;
            ctmp = ctag(mm2) ;
            dntmp = movmean(dntmp',Nrm) ;
            dtmp = movmean(dtmp',Nrm) ;
            ttmp = movmean(ttmp',Nrm) ;
            ctmp = movmean(ctmp',Nrm) ;

            dntmp_b{jj,ii} = dntmp;
            dtmp_b{jj,ii} = dtmp;
            ttmp_b{jj,ii} = ttmp;
            ctmp_b{jj,ii} = ctmp;
            
%             Dtag = [Dtag dtmp'] ;
%             Ttag = [Ttag ttmp'] ;
%             Ctag = [Ctag ctmp'] ;
%             DDrbr = [DDrbr interp1(dn,D,dntmp)'] ;
%             Trbr = [Trbr interp1(dn,T,dntmp)'] ;
%             Crbr = [Crbr interp1(dn,C,dntmp)'] ;
            
            if ii == 1
                figure(4)
                plot([1 1]*mean(dndn,'omitnan'),[0 50],'k') ;
            end

        
            Cmaxtag = max([Cmaxtag max(ctmp)]) ;
            
        end
        
    
    end
end
%%
ii=1

figure('Position',[444          92        1724        1225]), 
for il = 1:5
        for ic = 1:5 
            jj = 5*(il - 1) + ic ;
            subplot(5,5,jj)
            plot(tt_a{jj,ii},-dd_a{jj,ii},'--r','linewidth',1)
            hold on
            plot(ttmp_a{jj,ii},-dtmp_a{jj,ii},'--b','linewidth',1)
            plot(tt_b{jj,ii},-dd_b{jj,ii},'r','linewidth',1)
            plot(ttmp_b{jj,ii},-dtmp_b{jj,ii},'b','linewidth',1)
            if jj==21
                xlabel('Temp, degC')
                ylabel('Depth, m')
            end
            if il == 1
                set(gca,'xlim',[10 16])
            elseif il == 2
                set(gca,'xlim',[9 14])
            elseif il == 3
                set(gca,'xlim',[20 21])
            elseif il == 4
                set(gca,'xlim',[21.3 21.7])
            elseif il == 5
                set(gca,'xlim',[21.4 21.8])
            end

            if jj==1
                lg0 = legend('RBR down','Tag down','RBR up','Tag up');
                set(lg0,'position',[0.9090    0.8327    0.0721    0.0636])
            end
%             title([dcsv2(ii).name(1:6),',diff= ',num2str(mean(ttmp_a{jj,ii})-mean(tt_a{jj,ii}))])
            title([dcsv2(ii).name(1:6),', ',num2str(mean(ttmp_a{jj,ii})-mean(tt_a{jj,ii})),', ',...
                num2str(mean(ttmp_b{jj,ii})-mean(tt_b{jj,ii}))])
            set(gca,'fontsize',12)
        end
end
%%
figure('Position',[444          92        1724        1225]), 
for il = 1:5
        for ic = 1:5 ;
            jj = 5*(il - 1) + ic ;
            subplot(5,5,jj)
            plot(cc_a{jj,ii},-dd_a{jj,ii},'--r','linewidth',1)
            hold on
            plot(ctmp_a{jj,ii},-dtmp_a{jj,ii},'--b','linewidth',1)
            plot(cc_b{jj,ii},-dd_b{jj,ii},'r','linewidth',1)
            plot(ctmp_b{jj,ii},-dtmp_b{jj,ii},'b','linewidth',1)
            if jj==21
                xlabel('Conductivity, mS/cm')
                ylabel('Depth, m')
            end

            if il == 1
                set(gca,'xlim',[35 40])
            elseif il == 2
                set(gca,'xlim',[35 40])
            elseif il == 3
                set(gca,'xlim',[26 34])
            elseif il == 4
                set(gca,'xlim',[5 23])
            elseif il == 5
                set(gca,'xlim',[5.5 8])
            end

            if jj==1
                lg0 = legend('RBR down','Tag down','RBR up','Tag up');
                set(lg0,'position',[0.9090    0.8327    0.0721    0.0636])
            end
            title([dcsv2(ii).name(1:6),', ',num2str(mean(ctmp_a{jj,ii})-mean(cc_a{jj,ii})),', ',...
                num2str(mean(ctmp_b{jj,ii})-mean(cc_b{jj,ii}))])
            set(gca,'fontsize',14)
        end
end
%%
for il = 1:5
        for ic = 1:5 ;
            jj = 5*(il - 1) + ic ;
            dcc(jj) = mean(cc_a{jj})-mean(ctmp_a{jj});
        end
end
figure, plot(dcc)
%%
%     figure(2)
%     subplot(5,2,9)
%     xlabel('Temperature (^oC)')
%     subplot(5,2,10)
%     xlabel('Temperature (^oC)')
%     fnT = ['OrnitelaTests_YaquinaBay_202006_Temperature_' num2str(SN(ii)) '.jpg'] ;
%     eval(['print -djpeg90 -r300 ' fnT]) ;
%     figure(3)
%     subplot(5,2,9)
%     xlabel('Conductivity (mS/cm)')
%     subplot(5,2,10)
%     xlabel('Conductivity (mS/cm)')
%     fnC = ['OrnitelaTests_YaquinaBay_202006_Conductivity_' num2str(SN(ii)) '.jpg'] ;
%     eval(['print -djpeg90 -r300 ' fnC]) ;
% 
% end

figure(4)
datetick('x') ;
axis([dnmin dnmax 9 35])
title('Tag Interior Temperature vs. Time (Yaquina Bay, 06/25/2020)') ;
ylabel('T (^oC)') 
xlabel('Time') ;
print  -djpeg90 -r300 OrnitelaTests_YaquinaBay_202006_TagInteriorTemp.jpg

figure(5)
clf
% plot_size(1,5,15,5) ;
figure(6)
clf
% plot_size(1,5,15,5) ;
dn1 = datenum('26-Jun-2020 01:29:11') ;
dn2 = datenum('26-Jun-2020 01:38:30') ;
for ii = 1:4 ;
    nn = find((dn10Hz{ii}>=dn1).*(dn10Hz{ii}<=dn2)); 
    dntmp = dn10Hz{ii}(nn) ;
    Ctmp = C10Hz{ii}(nn) ;
    Ttmp = T10Hz{ii}(nn) ;
    figure(5)
    plot(dntmp,Ttmp,[co(ii) '.']) ;
    hold on
    figure(6)
    plot(dntmp,Ctmp,[co(ii) '.']) ;
    hold on
end

figure(5)
datetick('x') ;
axis([dn1 dn2 20.2 22.5])
title('Temperature vs Time (Yaquina Bay, Site 6, 06/25/2020)','fontsize',18) ;
ylabel('T (^oC)','fontsize',18) 
xlabel('Time','fontsize',18) ;
print  -djpeg90 -r300 OrnitelaTests_YaquinaBay_202006_TvsTime.jpg

figure(6)
datetick('x') ;
axis([dn1 dn2 4 16])
title('Conductivity vs Time (Yaquina Bay, Site 6, 06/25/2020)','fontsize',18) ;
ylabel('C (mS/cm)','fontsize',18) 
xlabel('Time','fontsize',18) ;
print  -djpeg90 -r300 OrnitelaTests_YaquinaBay_202006_CvsTime.jpg



