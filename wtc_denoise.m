% Script of the Wavelet-based method for removing fNIRS global physiological noise
% Please cite the paper: Duan et al., BOE, 2018, Wavelet-based method for removing global physiological noise in functional near-infrared spectroscopy 
% Input: data
% Output: denoised_data

raw_data = data;
threshold = 0.5;

N = size(raw_data,1);
N_ch = size(raw_data,2);

denoised_data = zeros(N, N_ch);
for cur_ch = 1:N_ch
    w_mat = zeros(101,N);
    for k = setdiff(1:N_ch, cur_ch)
        [wcoh,~,F,coi] = wcoherence(raw_data(:,cur_ch),raw_data(:,k),10,'VoicesPerOctave',10,'NumOctaves',10);
        w_binary = wcoh > 0.71;  % significance threshold. For speed, use 0.71 from Grinsted et al., 2004; Monte Calo can also be used, but slow. 
        w_mat = w_mat + w_binary;
    end
    
    mask = w_mat < N_ch * threshold;
    
    [cfs,f] = cwt(raw_data(:,cur_ch),'amor',10,'VoicesPerOctave',10,'NumOctaves',10);
    masked_cfs = cfs .* mask;
    xrec = icwt(masked_cfs);
    denoised_data(:,cur_ch) = xrec';
end