function X = applynoise_dep(A, p)
    d = length(A);
    trA = real(trace(A)); % for unitary U, tr UaU' = tr U'Ua = tr a
    X = (1-p) * A + (p/d)*trA*eye(d); 
end