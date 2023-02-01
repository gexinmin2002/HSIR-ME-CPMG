%%
%Function to add noise to echo trains
%%snr is the signal-to-noise level
function AtmTwnWithNoise=okaddNoise(AtmTwn,snr)
[row_AtmTwn,col_AtmTwn]=size(AtmTwn);
AtmTwnNoise=randn(row_AtmTwn,col_AtmTwn);
pSignal=sum(sum(AtmTwn.*AtmTwn))/(row_AtmTwn*col_AtmTwn);
pNoise=sqrt(pSignal/(10^(snr/10)));
AtmTwnWithNoise=AtmTwn+AtmTwnNoise*pNoise;