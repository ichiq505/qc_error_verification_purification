function CNOT = CNOT_ij_n(ctrl,target,num_qubit)
    % ctrl, target : ,1,2,3,...n
    
    X = [0,1;1,0];
    
    unit_ctrl0 = repmat([1,0], 2^(num_qubit-ctrl),1); vec_unit_ctrl0 = unit_ctrl0(:); clear unit_ctrl0
    vec_ctrl0 = repmat(vec_unit_ctrl0, 2^(ctrl-1),1); 
    CNOT = diag(vec_ctrl0); clear vec_ctrl0
    
    if ctrl < target
        diff_nt = num_qubit - target;
        X_rem = [zeros(2^diff_nt), eye(2^diff_nt); eye(2^diff_nt), zeros(2^diff_nt)];
        cX_rem = kron([0,0;0,1],blkdiagrep(X_rem, target-ctrl));
        CNOT_1 = blkdiagrep(cX_rem,ctrl);
        CNOT = CNOT + CNOT_1;
    else
        diff_nt = ctrl - target;
        X_mid = [zeros(2^diff_nt), eye(2^diff_nt); eye(2^diff_nt), zeros(2^diff_nt)];
        c_rem = repmat([0,1], 2^(num_qubit-ctrl),1); vec_c_rem = c_rem(:); clear c_rem
        c_rem = diag(vec_c_rem);
        CNOT_1 = blkdiagrep(kron(X_mid,c_rem),target);
        CNOT = CNOT + CNOT_1;
    end
end