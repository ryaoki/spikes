

function spikeStruct = loadKSdir(ksDir, varargin)
% modified to load params.sample_rate
% last updated @ 210601


if ~isempty(varargin)
    params = varargin{1};
else
    params = [];
end

if ~isfield(params, 'excludeNoise')
    params.excludeNoise = true;
end
if ~isfield(params, 'loadPCs')
    params.loadPCs = false;
end

% load spike data

spikeStruct = loadParamsPy(fullfile(ksDir, 'params.py'));

if isfield(params, 'sample_rate')
    spikeStruct.sample_rate = params.sample_rate; % overwrite default acquisition rate by empirical acquisition rate
end

ss = readNPY(fullfile(ksDir, 'spike_times.npy'));
st = double(ss)/spikeStruct.sample_rate;
spikeTemplates = readNPY(fullfile(ksDir, 'spike_templates.npy')); % note: zero-indexed

if exist(fullfile(ksDir, 'spike_clusters.npy'))
    clu = readNPY(fullfile(ksDir, 'spike_clusters.npy'));
else
    clu = spikeTemplates;
end

tempScalingAmps = readNPY(fullfile(ksDir, 'amplitudes.npy'));

if params.loadPCs
    pcFeat = readNPY(fullfile(ksDir,'pc_features.npy')); % nSpikes x nFeatures x nLocalChannels
    pcFeatInd = readNPY(fullfile(ksDir,'pc_feature_ind.npy')); % nTemplates x nLocalChannels
else
    pcFeat = [];
    pcFeatInd = [];
end

cgsFile = '';
if exist(fullfile(ksDir, 'cluster_groups.csv')) 
    cgsFile = fullfile(ksDir, 'cluster_groups.csv');
end
if exist(fullfile(ksDir, 'cluster_group.tsv')) 
   cgsFile = fullfile(ksDir, 'cluster_group.tsv');
end 
if ~isempty(cgsFile)
    [cids, cgs] = readClusterGroupsCSV(cgsFile);

    if params.excludeNoise
        noiseClusters = cids(cgs==0);

        st = st(~ismember(clu, noiseClusters));
        spikeTemplates = spikeTemplates(~ismember(clu, noiseClusters));
        tempScalingAmps = tempScalingAmps(~ismember(clu, noiseClusters));        
        
        if params.loadPCs
            pcFeat = pcFeat(~ismember(clu, noiseClusters), :,:);
            %pcFeatInd = pcFeatInd(~ismember(cids, noiseClusters),:);
        end
        
        clu = clu(~ismember(clu, noiseClusters));
        cgs = cgs(~ismember(cids, noiseClusters));
        cids = cids(~ismember(cids, noiseClusters));
        
        
    end
    
else % skip phy
    clu = spikeTemplates;
    
    cids = unique(spikeTemplates);
    cgs = 3*ones(size(cids));
end
    

coords = readNPY(fullfile(ksDir, 'channel_positions.npy'));
ycoords = coords(:,2); xcoords = coords(:,1);
temps = readNPY(fullfile(ksDir, 'templates.npy'));

winv = readNPY(fullfile(ksDir, 'whitening_mat_inv.npy'));

% idx for clusters to template: this is necessary if you merge/split any cluster in phy
for iUni = 1: length(cids)
    cluToTemp(iUni) = mode(spikeTemplates(clu==cids(iUni)));
end




spikeStruct.st = st;
spikeStruct.spikeTemplates = spikeTemplates;
spikeStruct.clu = clu;
spikeStruct.tempScalingAmps = tempScalingAmps;
spikeStruct.cgs = cgs;
spikeStruct.cids = cids;
spikeStruct.xcoords = xcoords;
spikeStruct.ycoords = ycoords;
spikeStruct.temps = temps;
spikeStruct.winv = winv;
spikeStruct.pcFeat = pcFeat;
spikeStruct.pcFeatInd = pcFeatInd;
spikeStruct.cluToTemp = cluToTemp;