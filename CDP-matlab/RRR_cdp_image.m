function RRR_cdp_image(b, y, sol, D)

    %% input parameters
    beta = 0.5 ;
    epsilon = 1.e-5 ;
    iter_max = 100000 ;
    L = size(D, 1) ;
    height = size(D, 2) ;
    width = size(D, 3) ;

    %% The inverse of the diagonal terms of 1/n * Hermitian(A) * A
    inv_diag = zeros(height, width) ;
    for i = 1:L
        for j = 1:height
            for k = 1:width
                inv_diag(j, k) = inv_diag(j, k) + abs(D(i, j, k))^2 ;
            end
        end
    end
    inv_diag = 1 ./ inv_diag ;

    iter_arr = [] ;
    err_arr = [] ;
    cur_sol = zeros(height, width, 3) ;
    for iter = 1:iter_max
        rms_diff = 0. ;
        norm_y = 0. ;
        for c = 1:3
            b_c = squeeze(b(:, :, :, c)) ;
            y_c = squeeze(y(:, :, :, c)) ;
            p1 = proj1(D, inv_diag, y_c, L, height, width) ;
            p2 = proj2(b_c, 2*p1 - y_c) ;
            delta_y = p2 - p1 ;
            for k = 1:L
                rms_diff = rms_diff + norm(squeeze(delta_y(k, :, :)))^2 ;
                norm_y = norm_y + norm(squeeze(y_c(k, :, :)))^2 ;
            end
            y_c = y_c + beta*delta_y ;
            y(:, :, :, c) = y_c ;
        end

        rms_diff = beta*sqrt(rms_diff / norm_y) ;
        if (iter < 10 || mod(iter, 10) == 0)
            norm1 = 0. ;
            norm2 = 0. ;
            for c = 1:3
                x_c = reshape(squeeze(sol(:, :, c)), height*width, 1) ;
                y_c = squeeze(y(:, :, :, c)) ;
                p1 = proj1(D, inv_diag, y_c, L, height, width) ;
                cur_sol(:, :, c) = Ainv(D, inv_diag, p1, L, height, width) ;
                x_cur = reshape(squeeze(cur_sol(:, :, c)), height*width, 1) ;
                norm1 = norm1 + norm(x_c - x_cur*sign(x_cur'*x_c))^2 ;
                norm2 = norm2 + norm(x_c)^2 ;
            end
            filename = ['reconst-cornell-iter', int2str(iter), '.jpg'] ;
            cur_img = uint8(round(abs(cur_sol))) ;
            imwrite(cur_img, filename) ;
            err_val = sqrt(norm1/norm2) ;
            iter_arr = [iter_arr, iter] ;
            err_arr = [err_arr, err_val] ;
            [err_val, rms_diff]
            if (rms_diff < epsilon)
                break ;
            end
        end
    end

    fp = fopen('phasing-error.dat', 'w') ;
    for i = 1:length(err_arr)
        fprintf(fp, '%d %1.5e\n', iter_arr(i), err_arr(i)) ;
    end
    fclose(fp) ;
end

%% pseudo-inverse of A
function x1 = Ainv(D, inv_diag, y, L, height, width)
    x1 = zeros(height, width) ;
    for l = 1:L
        x1 = x1 + squeeze(conj(D(l, :, :))) .* ifft2(squeeze(y(l, :, :))) ;
    end
    x1 = x1 .* inv_diag ;
end

%% projections
function y1 = proj1(D, inv_diag, y, L, height, width)
    x1 = Ainv(D, inv_diag, y, L, height, width) ;
    y1 = zeros(L, height, width) ;
    for i = 1:L
        y1(i, :, :) = fft2(squeeze(D(i, :, :)) .* x1) ;
    end
end

function y2 = proj2(b, y)
    y2 = b .* sign(y) ;
end
