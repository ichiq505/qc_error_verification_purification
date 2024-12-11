function [purified_fids, probs_succ] = state_purify_py(fid, error_meas, error_cx, max_add)
    state_0 = [1,0;0,0];
    % meas_0 = [1,0;0,0];
    % meas_1 = [0,0;0,1];
    
    dep_st = fid2dep(fid);
    
    dep_meas = 2*error_meas;

    noisy_state = applynoise_dep(state_0, dep_st);
    
    num_additional = 0:max_add;
    
    purified_fids = zeros(size(num_additional));
    probs_succ = zeros(size(num_additional));
    for num_a=1:length(num_additional)
        [purified_state, prob_succ] = state_purify(noisy_state,num_additional(num_a),dep_st,error_cx,dep_meas);
        purified_fids(num_a) = purified_state(1);
        probs_succ(num_a) = prob_succ;
    end

    % subplot(211), hold on
    % plot(num_additional,purified_fid,'o-')
    % subplot(212), hold on
    % plot(num_additional,probs_succ,'o-')
    % cx_labels{error_cx} = "$\epsilon=$"+string(errors_cx(error_cx));
    % 
    % subplot(211)
    % legend(cx_labels,'Interpreter','latex','Location','southeast'), ylabel('$f^{(n)}$','Interpreter','latex'), axis([0 N 0.9 1]), grid on, title(["state preparation error $(1-f)=$"+string(1-fid),"measurement error $q=$"+string(q)],'Interpreter','latex') %,"CNOT gate error $\epsilon=$"+string(error_cx)
    % subplot(212)
    % legend(cx_labels,'Interpreter','latex','Location','northeast'), ylabel('$p_{succ}^{(n)}$','Interpreter','latex'), axis([0 N -inf 1]), grid on
end