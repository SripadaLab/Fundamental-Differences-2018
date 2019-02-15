function [latent, latentLow, latentHigh, latentShuffle] = perm_pca_timeseries(SubjDir,TCTemplate,RunDir, nShuffle)

if ~exist('nShuffle', 'var')
        nShuffle = 100;
end

%under new permutation scheme we need to pass in subject IDs and path
%template for ROI timecourses.
n = size(SubjDir,1);
t = 4800*.9; %minimum length to trim to due to different censored lengths
r = 264; %power nodes

ts = zeros(n,t,r);
dur = [];
for iSubject = 1:n
    tcf = strrep(TCTemplate,'[Subject]',SubjDir(iSubject,:));
    tc = [];
    for iRun = 1:size(RunDir,1)
        tcfr = strrep(tcf,'[Run]',RunDir{iRun});
        temp = load(tcfr);
        tc = [tc;temp.roiTC];
        fprintf(1,'%d',size(temp.roiTC,1));
        dur = [dur;size(temp.roiTC,1)];
    end
    nt = size(tc,1);
    idx = sort(randsample(nt,t,0));
    ts(iSubject,:,:) = tc(idx,:);
end

%baseline
x = zeros(n,(r*(r-1))/2);
for iSubject = 1:n
    c = corr(squeeze(ts(iSubject,:,:)));
    x(iSubject,:) = mc_flatten_upper_triangle(c);
end
[~,~,d] = pca(x);
latent = sort(d,'descend');

xShuffle = x;
latentShuffle = zeros(length(latent), nShuffle);
parfor iShuffle = 1:nShuffle
    xShuffle = zeros(n,(r*(r-1))/2);
    perms = zeros(n,r);
    for iR = 1:r
        perms(:,iR) = randperm(n);
    end
    
    for iSubject = 1:n
        tempts = zeros(t,r);
        for iR = 1:r
            tempts(:,iR) = ts(perms(iSubject,iR),:,iR);
        end
        c = corr(tempts);
        xShuffle(iSubject,:) = mc_flatten_upper_triangle(c);
    end
    
    [~,~,d] = pca(xShuffle);
    latentShuffle(:,iShuffle) = sort(d,'descend')';
end

latentHigh = quantile(latentShuffle', 1-0.025);
latentLow = quantile(latentShuffle', 0.025);
