function Y = ToepPhi(c,X)
N = size(X,2);
xhat = fft2(c) .* fft2(X,2*N,2*N);
xmat = ifft2(xhat);
Y = xmat(1:N,1:N);
end
