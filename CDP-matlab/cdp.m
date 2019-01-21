L = 3 ;
stanford = double(imread('stanford.jpg')) ;
cornell = double(imread('cornell.jpg')) ;
height = size(cornell, 1) ;
width = size(cornell, 2) ;

%% generate random coded diffraction pattern
D = zeros(L, height, width) ;
for i = 1:L
    D(i, :, :) = exp(1i*floor(rand(height, width)*4)*pi/2) ;
end

b = zeros(L, height, width, 3) ;    % data
y = zeros(L, height, width, 3) ;    % initial guess
for c = 1:3
    for i = 1:L
        phase_mask = squeeze(D(i, :, :)) ;
        b(i, :, :, c) = abs(fft2(phase_mask .* cornell(:, :, c))) ;
        y(i, :, :, c) = fft2(phase_mask .* stanford(:, :, c)) ;
    end
end

RRR_cdp_image(b, y, cornell, D) ;
