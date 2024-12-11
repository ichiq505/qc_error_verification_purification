%% [purified_state, prob_succ] = state_purify(noisy_state, num_additional, error_additional, error_cx, error_meas)
%     noisy_state: state (density matrix) to be purified
%     num_additional: number of additional qubits for purification
%     error_additional: (list) state preparation error rates of additional qubits
%                * if length(error_additional) < num_additional, error rate for all additional qubits 
%                  is set to error_additional(1).
%                * depolarizing noise: dep = 2 - 2f
%     error_cx: CNOT error rate (depolarizing noise)
%     error_meas: measurement error rate (depolarizing noise: dep = 2q)
%
function [purified_state, prob_succ] = state_purify(noisy_state, num_additional, error_additional, error_cx, error_meas)
%     nargin
    switch nargin
        case 1
            disp("'num_additional' is not given. No purification.")
            purified_state = noisy_state;
            prob_succ = 1;
            return
        case 2
            disp('num_additional : '+string(num_additional)+' / no error list.')
%             error_cx = 0;
%             error_meas = 0;
        case 3
            disp('num_additional : '+string(num_additional)+' / error list for additional qubits is given.')
%             error_cx = 0;
%             error_meas = 0;
        case 4
            disp('num_additional : '+string(num_additional)+' / CNOT error is given.')
%             error_meas = 0;
        case 5
            disp('num_additional : '+string(num_additional)+' / measurement error is given.')
        otherwise
            error("Inputs are 'state','number of additional qubits', 'error rate for additional state preparations', 'error rate for CNOT gate'");
    end
        
    if ~(num_additional>0)
        disp('==>> num_additional is 0. No purification.')
        purified_state = noisy_state;
        prob_succ = 1;
        return
    else
        % varargin check: if yes, different error rate for additional qubits
        num_errors = length(error_additional);
        if num_errors==0
            disp('==>> the list of errors is empty. the errors for all additional qubits is set to 0.')
            error_add = zeros(num_additional,1);
        elseif num_errors == num_additional
            disp('==>> size of list of errors :'+string(num_errors))
            error_add = error_additional;
        else
            disp('==>> size of list of errors is not equal to '+string(num_errors))
            disp('==>> error rate for all additional qubits is set to error_additional(1)')
            error_add = repmat(error_additional(1),num_additional,1);
        end
    end

    state_0 = [1,0;0,0];
    dim_sys = 2; %length(noisy_state);
    dim_add = 2^num_additional;
    noisy_add_state = eye(dim_add);
    for n=1:num_additional
%         noisy_add_state_k = @(error) (1-error)*state_0 + error*eye(2)/2;
        noisy_add_state(1:2^n,1:2^n) = kron(noisy_add_state(1:2^(n-1),1:2^(n-1)),applynoise_dep(state_0,error_add(n)));
    end
    
    CX_2 = eye(4); CX_2(3:4,:) = CX_2(4:-1:3,:); %(3)
    sub_state = zeros(dim_sys^2); %(3)
    whole_state = kron(noisy_state,noisy_add_state);
%     CX_n = eye(dim_sys*dim_add); %(1)
    for n=1:num_additional
%         CX_n = CNOT_ij_n(1,1+n,1+num_additional) * CX_n; %(1)

%         CX_1n = CNOT_ij_n(1,1+n,1+num_additional); %(2)
%         whole_state = applynoise_dep(CX_1n*whole_state*CX_1n',error_cx); %(2)
        
        for row_sys=0:1  %(3)
            for row_add=0:1  %(3)
                row_sub = row_sys*2+row_add+1; %(3)
                row_1n = row_sys*2^(num_additional+1-1)+row_add*2^(num_additional+1-(n+1))+1; %(3)
                for col_sys=0:1  %(3)
                    for col_add=0:1  %(3)
                        col_sub = col_sys*2+col_add+1; %(3)
                        col_1n = col_sys*2^(num_additional+1-1)+col_add*2^(num_additional+1-(n+1))+1; %(3)
                        sub_state(row_sub,col_sub) = whole_state(row_1n,col_1n); %(3)
                    end
                end
            end
        end
        sub_state = applynoise_dep(CX_2*sub_state*CX_2',error_cx); %(3)
        for row_sys=0:1  %(3)
            for row_add=0:1  %(3)
                row_sub = row_sys*2+row_add+1; %(3)
                row_1n = row_sys*2^(num_additional+1-1)+row_add*2^(num_additional+1-(n+1))+1; %(3)
                for col_sys=0:1  %(3)
                    for col_add=0:1  %(3)
                        col_sub = col_sys*2+col_add+1; %(3)
                        col_1n = col_sys*2^(num_additional+1-1)+col_add*2^(num_additional+1-(n+1))+1; %(3)
                        whole_state(row_1n,col_1n) = sub_state(row_sub,col_sub); %(3)
                    end
                end
            end
        end
    end
%     whole_state = applynoise_dep(CX_n * whole_state * CX_n',error_cx); %(1)
    
    if error_meas ~= 0
        q = error_meas / 2; % dep -> q
        coefs = zeros(dim_add,1);
        for i=1:2^num_additional 
            num_ones = nnz(double(dec2bin(i-1))-double('0'));
            coefs(i) = q^num_ones * (1-q)^(num_additional-num_ones);
        end
        meas_add = diag(coefs);
    else % error_meas = 0
        meas_add = diag(repmat(([1; zeros(dim_add-1,1)]),2,1));
    end
    meas_all = kron(eye(dim_sys),meas_add);
    purified_state = PartialTrace(whole_state*meas_all,2,[dim_sys,dim_add]);
    prob_succ = real(trace(purified_state));
    purified_state = purified_state / prob_succ;
    disp(' ')
end

%% Helper functions
