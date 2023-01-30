function [c,A] = getc(N,h,getA)

for jj = 1:N
    for kk = 1:N
        absj = sqrt( (jj-1)^2 + (kk-1)^2 );
        if jj == 1 && kk==1
            VecColumn(1) = 0;
        else
            VecColumn(N*(jj-1) + kk,1) = 1i/4*besselh(0,absj*h);
        end
    end
end
for ii = 1:N
    for jj = 1 : N
        absj = sqrt( (ii-1)^2 + (jj-1)^2 );
        if ii == 1 && jj==1
            ColumnForT(ii,jj) = 0;
        else
            ColumnForT(ii,jj) = 1i/4*besselh(0,absj*h);
        end
    end
end
for jj = 1:N
    for kk = 1:N
        absj = sqrt( (jj-1)^2 + (kk-1)^2 );
        if jj == 1 && kk==1
            VecRow(1) = 0;
        else
            VecRow(1,N*(jj-1) + kk) = 1i/4*besselh(0,absj*h);
        end
    end
end

if nargin == 3 && getA == 1
    % create first row blocks
    RowMats = {};
    for jj = 1:N
        RowMats{end+1} = toeplitz(ColumnForT(jj,:),[ColumnForT(jj,1), VecRow(N*(jj-1)+2:N*jj)]);
    end
    n=length(RowMats);
    A = cell2mat(RowMats(toeplitz(1:n)));
else
    A = [];
end
c = zeros(2*N,2*N);
c(1:N,1:N) = ColumnForT;
for jj = 1:N
    c(N+2:end,jj) = fliplr(VecRow(N*(jj-1)+2:N*jj));
end
c(:,N+2:end) = fliplr(c(:,2:N));
end


