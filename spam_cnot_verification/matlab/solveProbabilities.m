function solution = solveProbabilities(Prob00, Prob01, Prob10, Prob11)
    % Define the function that returns the array of equations
    equations = @(v) [
        (1-v(3))*((1-v(1))^2*(1-v(2))^2 + (1-v(1))*v(1)*(1-v(2))*v(2) + v(1)^2*v(2)*(1-v(2)) + v(1)*(1-v(1))*v(2)^2) + v(3)/4 - Prob00;
        (1-v(3))*((1-v(1))^2*(1-v(2))*v(2) + (1-v(1))*v(1)*(1-v(2))^2 + v(1)^2*v(2)^2 + v(1)*(1-v(1))*v(2)*(1-v(2))) + v(3)/4 - Prob01;
        (1-v(3))*((1-v(1))^2*v(2)*(1-v(2)) + (1-v(1))*v(1)*v(2)^2 + v(1)^2*(1-v(2))^2 + v(1)*(1-v(1))*(1-v(2))*v(2)) + v(3)/4 - Prob10;
        (1-v(3))*((1-v(1))^2*v(2)^2 + (1-v(1))*v(1)*v(2)*(1-v(2)) + v(1)^2*(1-v(2))*v(2) + v(1)*(1-v(1))*(1-v(2))^2) + v(3)/4 - Prob11
    ];

    % Initial guesses for p and q
    initial_guess = [0, 0, 0];

    % Use fsolve to solve the equations
    options = optimoptions('fsolve', 'Display', 'none'); % Set Display to 'none'
    % [solution, ~, exitflag, ~] = fsolve(equations, initial_guess, options);
    solution = fsolve(equations, initial_guess, options);

    disp('Solution found:');
    disp(solution);

    % Check if solution was successful
    % if exitflag ~= 1
    %     disp('The solver did not converge to a solution.');
    % else
    %     disp('Solution found:');
    %     disp(solution);
    % end
end