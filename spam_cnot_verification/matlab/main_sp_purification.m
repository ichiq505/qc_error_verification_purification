clear, clf, clc

state_0 = [1,0;0,0];
meas_0 = [1,0;0,0];
meas_1 = [0,0;0,1];

fid = 0.95;
dep = fid2dep(fid);

errors_cx = [0 .01 .03 .05 .1];
error_meas = dep; 
q = error_meas / 2; % working...

noisy_state = applynoise_dep(state_0, dep);

N = 5;
num_additional = 0:N;

cx_labels = cell(length(errors_cx),1);
for error_cx=1:length(errors_cx)
    purified_fid = zeros(size(num_additional));
    probs_succ = zeros(size(num_additional));
    for num_a=1:length(num_additional)
        [purified_state, prob_succ] = state_purify(noisy_state,num_additional(num_a),dep,errors_cx(error_cx),dep);
        purified_fid(num_a) = purified_state(1);
        probs_succ(num_a) = prob_succ;
    end
    subplot(211), hold on
    plot(num_additional,purified_fid,'o-')
    subplot(212), hold on
    plot(num_additional,probs_succ,'o-')
    cx_labels{error_cx} = "$\epsilon=$"+string(errors_cx(error_cx));
end

subplot(211)
legend(cx_labels,'Interpreter','latex','Location','southeast'), ylabel('$f^{(n)}$','Interpreter','latex'), axis([0 N 0.9 1]), grid on, title(["state preparation error $(1-f)=$"+string(1-fid),"measurement error $q=$"+string(q)],'Interpreter','latex') %,"CNOT gate error $\epsilon=$"+string(error_cx)
subplot(212)
legend(cx_labels,'Interpreter','latex','Location','northeast'), ylabel('$p_{succ}^{(n)}$','Interpreter','latex'), axis([0 N -inf 1]), grid on


%% Helper functions
