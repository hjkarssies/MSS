function [x, active_set] = incaActiveSet(Q, c, A, b, x0)
% incaActiveSet - Active Set QP Solver
%
% Solves: min 0.5 * x' Q x + c' x
%         s.t. A x <= b
%
% Inputs:
%   Q, c   - Cost function terms
%   A, b   - Inequality constraint matrix and bounds
%   x0     - Initial guess
%
% Output:
%   x            - Optimal solution
%   active_set   - Final set of active constraints

max_iter = 20;
tol = 1e-6;

x = x0;
active_set = find(abs(A * x - b) < tol);  % Start with tight constraints

for k = 1:max_iter
    if isempty(active_set)
        dx = -Q \ c;
    else
        Aeq = A(active_set, :);
        beq = b(active_set);

        try
            % Solve equality-constrained QP
            KKT = [Q, Aeq'; Aeq, zeros(length(active_set))];
            rhs = [-c; beq - Aeq * x];
            sol = KKT \ rhs;
            dx = sol(1:length(x));
            lambda = sol(length(x)+1:end);
        catch
            warning('KKT system is singular at iteration %d', k);
            break;
        end
    end

    x_new = x + dx;

    % Find most violated inequality
    viol = find(A * x_new - b > tol);
    if ~isempty(viol)
        d = dx;
        alpha = 1;
        for j = viol'
            a = A(j,:);
            if a * d > 0
                alpha_j = (b(j) - a * x) / (a * d);
                if alpha_j < alpha
                    alpha = alpha_j;
                    worst = j;
                end
            end
        end
        x = x + alpha * dx;
        if ~ismember(worst, active_set)
            active_set = [active_set; worst];
        end
    else
        if isempty(active_set)
            break;
        end
        if all(lambda >= -tol)
            x = x_new;
            break;
        else
            [~, idx] = min(lambda);
            active_set(idx) = [];
            x = x_new;
        end
    end
end

end
