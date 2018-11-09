clear all
%close all
folder = uigetdir();
set(0,'defaultAxesFontName','Arial');
set(0,'defaultAxesFontSize',20);
allL = [];
allF = [];
allEx = [];
ColOrd = get(gca,'ColorOrder');
leading_edge=@(x,v)((x(:) > v & circshift(x(:),[1 0]) < v));
%%
d = dir([folder filesep '*-motion.csv']);
d = d(3);
for i = 1:length(d)
    ex = strrep(d(i).name, '-motion.csv','');
    disp(ex);
    motion = dlmread([folder filesep d(i).name],',',1,0);
    t = motion(:,1);
    y = motion(:,3);
    % plot(t,y);pause;

    % load the oscilloscope trace
    d2 = dir([ folder filesep 'psdata' filesep ex filesep '*.mat']);
    if isempty(d2)
        error('Folder ''%s'' is empty', [folder filesep 'psdata' filesep ex filesep '*.mat']);
    end
    fire = [];
    pulse = [];
    for j = 1:length(d2)
        psname = [folder filesep 'psdata' filesep ex filesep d2(j).name];
        psdata = load(psname);
        fire = [fire; psdata.A];
        pulse = [pulse; psdata.B];
    end
    time = psdata.Tinterval * (1:length(fire))';

    % compare time stamps and fire signal
    % fire signal
    sig1 = leading_edge(fire,2.5);
    % timestamps
    sig2 = zeros(size(time));
    sig2(round(t(t<max(time))/psdata.Tinterval)) = 1;
    %sig2 = sig2(1:numel(sig1));
    figure(10); subplot(numel(d),1,i)
    plot(time,cumsum(leading_edge(sig1,0.5))-cumsum(sig2));
    title('discrepency between fire and timestamps');
    ylabel('number of frame')
    xlabel('time [s]')
    %
    % detect leading edge of both ttl signals
    tfire = time(leading_edge(-fire,-2.5));
    tpulse = time(leading_edge(pulse,2.5));

    % correct missing TTL signals
    % Some TTL signal are missing due to the gap between segments recored
    % by the picoscope
    % Guess the correct time intercal using the difference between the two
    % first pulses.
    dt = tpulse(2) - tpulse(1);
    pulseFreq = 1/dt;
    fprintf('Pulse freq %f Hz\n', pulseFreq);
    tpulsec = tpulse(1) + (0:round(max(t)/dt))*dt;

    % compute the acquisition frequency from the ttl signal
    T = tfire(end) - tfire(1);
    dt = median(diff(tfire));
    fprintf('Fire freq %f Hz\n', 1/dt);
    tfirec = tfire(1) + (0:round(T/dt))*dt;

    fprintf('Timestamp freq %f Hz (%.2fms)\n', 1/mean(diff(t)), 1000*mean(diff(t)));


    % synchronize the oscilloscope trace and the nd2 timestamps by shifting
    % the time stamps of the movie to the first fire signal from the camera.
    % we use the time stamps from the nd2 file as the TTL signal has too
    % many gaps in it to be reliable
    exposuretime = median(time(leading_edge(-fire,-2.5)) - time(leading_edge(fire, 2.5)));
    fprintf('Exposure time (fire) %.2fms\n', 1000*exposuretime);
    t = t + tfire(1) - t(1) + dt;

    % motion detection
    motionthreshold = 2 * std(y) + median(y);   
    tmotion  = t(leading_edge(y, motionthreshold));
    
    figure(11); subplot(numel(d),1,i);
    plot(tpulse,ones(size(tpulse)),'.','markersize',20); hold on;
    plot(tpulsec,ones(size(tpulsec))+0.001,'.','markersize',20);
    plot(tfire,ones(size(tfire))+0.002,'.','markersize',20);
    plot(t,ones(size(t))+0.003,'.','markersize',20);
    plot(tfirec,ones(size(tfirec))+0.004,'.','markersize',20);
    plot(tmotion,ones(size(tmotion))+0.005,'.','markersize',20);   
    hold off; ylim([0.9 1.2])
    title('fire signal');
    legend('pulse','pulse corrected','fire recorded','fire timestamps','fire corrected','motion')
    ylabel('TTL')
    xlabel('Time [s]')



    % for each stimulation find the closest motion activity in the future
    % and deduce the latency L
    L = -ones(1,length(tpulsec)-1);

    %% Compute average curve
    figure(1);
    subplot(2*length(d),1,2*i-1);
    %subplot(2,2,i);
    % accumulators for storing the average curves
    acc0=[];
    acc1=[];
    acc2=[];
    cla;
    
    Dt = tpulsec(2) - tpulsec(1);
    for k = 1:length(tpulsec)
        % index for t and y in between two consecutives pulses
        tidx = find(t<=tpulsec(k),1,'last'):find(t<tpulsec(k)+Dt,1,'last');
        if ~isempty(tidx)
            l = find(y(tidx) > motionthreshold,1); % index of the next motion activity
            if ~isempty(l)
                L(k) = t(tidx(l)) - tpulsec(k); % latency between the two peaks
                plot(t(tidx) - tpulsec(k),y(tidx),'Color',[.7 .7 .7]);
                hold on;
                if isempty(acc0)
                    acc0 = zeros(length(tpulsec), 2*length(tidx));
                    acc1 = zeros(length(tpulsec), 2*length(tidx));
                    acc2 = zeros(length(tpulsec), 2*length(tidx));
                end
                acc0(k,1:length(tidx)) = 1;
                acc1(k,1:length(tidx)) = y(tidx);
                acc2(k,1:length(tidx)) = t(tidx) - tpulsec(k);
            end
        end
    end
    plot(sum(acc2) ./ sum(acc0), sum(acc1)./sum(acc0), 'k', 'linewidth', 4);        
    hold off
    axis([0 1 0 max(y)]);
    %title(sprintf('Stimulation frequency: %.2fHz Latency:%.2fms', pulseFreq, 1000*median(L(L~=-1))),'interpreter','none')
    xlabel('Latency [s]');
    ylabel('Motion [a.u.]');
    xlim([0 1])
    box off
    
    

    subplot(2*length(d),1,2*i);   
    datL = L(L~=-1);
    boxplot(datL,'orientation','horizontal')
    xlim([0 1])
    box off

    %% plot the overall signal
    idx = L > 0;
    for figid = [2,4]
        figure(figid);
        if figid==2
            %set(gcf,'Position',[0 0 1920 1080]);
            %subplot(length(d),1,i);
            subplot(numel(d),1,i)
        else
            %set(gcf,'Position',[0 0 1920 1080]);
            figure(figid);
            clf;
        end
        plot(t, y, '.-', 'markersize',5);
        hold on;
        for k = 1:length(idx)
            if idx(k) > 0
                plot([tpulsec(k) tpulsec(k)+L(k)],[-motionthreshold/10,-motionthreshold/10],'Color',ColOrd(3,:),'linewidth',4);
            end
        end
        plot(tpulsec,-motionthreshold/10*ones(size(tpulsec)),'.','Color',ColOrd(4,:),'MarkerSize',20);
        plot(tpulsec(idx)+L(idx),-motionthreshold/10*ones(size(L(idx))),'.','Color',ColOrd(2,:),'MarkerSize',20);


        hold off
        axis tight
        title(sprintf('Stimulation frequency: %.2fHz', pulseFreq),'interpreter','none')
        %title(ex,'interpreter','none')
        xlabel('Time [s]');
        ylabel('Motion');
        if figid==4
            figure(figid);
            print(sprintf('%s%straces-%d.eps',folder,filesep,i),'-depsc2')
        end
    end

    %legend('frame difference','pulse ttl','pulse corrected')
    fprintf(1,'latency median:%f ms std:%f ms n:%d\n', 1000*median(L(L~=-1)), 1000*std(L(L~=-1)), length(find(L~=-1) ));

    allL = [allL L(L~=-1)];
    allF = [allF pulseFreq * ones(size(L(L~=-1)))];
    allEx = [allEx i*ones(size(L(L~=-1)))];
end
%%
data = table(allEx',allF', allL','VariableNames',{'Sample','Freq','Latency'});
writetable(data,[folder filesep 'LatencyVsFreq.csv']);
figure(3);clf;
boxplot(allL,round(allF,1));
ylim([0 1.1*max(allL)]);
xlabel('Stimualtion frequency [Hz]');
ylabel('Latency [s]');
figure(5);
boxplot(allL,allEx)
xlabel('Sample');
ylabel('Latency [s]');

%%
figure(1)
print('-depsc2', [folder filesep 'average traces.eps']);
figure(2)
print('-depsc2', [folder filesep 'motion and stimualtion traces.eps']);
figure(3)
print('-depsc2', [folder filesep 'latency boxplot.eps']);
