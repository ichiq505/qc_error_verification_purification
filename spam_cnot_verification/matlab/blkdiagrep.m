function blkdiag_mtx = blkdiagrep(A,n)
    if n==0
        blkdiag_mtx = [];
        return
    end

    size_A = size(A);
    blkdiag_mtx = zeros(size_A * 2^(n-1));
    
    for k=1:2^(n-1)
        blkdiag_mtx(size_A(1)*(k-1)+1:size_A(1)*k,size_A(2)*(k-1)+1:size_A(2)*k) = A;
    end
end