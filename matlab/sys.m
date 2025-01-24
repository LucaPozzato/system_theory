clc

syms x1 x2 R L C a b 

% assumptions for symbolical analysis
assume([L == 11, C == 11, a == 0.04, b == 0.18]);

% system
dx1 = -(R*x1)/L+x2/L; % dx1/dt
dx2 = -x1/C+(a*x2)/C-(b*x2^3)/C; % dx2/dt

% solve for equilibria
[eq_x1, eq_x2] = solve([dx1 == 0, dx2 == 0], [x1, x2]);
equilibria = [eq_x1, eq_x2];
disp('Equilibria:');
disp(equilibria);

% evaluating stability with linearization method
% jacobian matrix
J = jacobian([dx1, dx2], [x1, x2]);
disp('Jacobian Matrix:');
disp(J);

% trace
tr_J = trace(J);
disp('Trace of Jacobian:');
disp(tr_J)

% determinant
det_J = det(J);
disp('Determinant of Jacobian')
disp(det_J)

R_tr_neg_eq = [sym(zeros(1, 10)); sym(zeros(1, 10))];

R_det_pos_eq = [sym(zeros(1, 10)); sym(zeros(1, 10))];

% evaluating stability in terms of R
for i = 1:length(eq_x1)
    % substitute the equilibrium point into the Jacobian
    J_eq = subs(J, {x1, x2}, {eq_x1(i), eq_x2(i)});
    
    % calculate the trace and determinant
    tr_J_eq = trace(J_eq);
    det_J_eq = det(J_eq);

    R_tr_neg = solve(tr_J_eq < 0, R, 'ReturnConditions', true, 'Real', true);
    R_tr_neg = subs(R_tr_neg.conditions, 'x', 'R');
    if length(R_tr_neg) > 1
        for k = 1:length(R_tr_neg)
            R_tr_neg_eq(i, k) = R_tr_neg(k);
        end
    else
        R_tr_neg_eq(i, 1) = R_tr_neg;
    end


    R_det_pos = solve(det_J_eq > 0, R, 'ReturnConditions', true, 'Real', true);
    R_det_pos = subs(R_det_pos.conditions, 'x', 'R');
    if length(R_det_pos) > 1
        for k = 1:length(R_det_pos)
            R_det_pos_eq(i, k) = R_det_pos(k);
        end
    else
        R_det_pos_eq(i, 1) = R_det_pos;
    end
end

eq_x1_existence = [sym(zeros(1, 10)); sym(zeros(1, 10))];

eq_x2_existence = [sym(zeros(1, 10)); sym(zeros(1, 10))];

for i = 1:length(eq_x1)
    x1_conditions = solve(imag(eq_x1(i)) == 0, R, 'ReturnConditions', true, 'Real', true);
    x1_conditions = subs(x1_conditions.conditions, 'x', 'R');
    if length(x1_conditions) > 1
        for k = 1:length(x1_conditions)
            eq_x1_existence(i, k) = x1_conditions(k);
        end
    else
        eq_x1_existence(i, 1) = x1_conditions;
    end

    x2_conditions = solve(imag(eq_x2(i)) == 0, R, 'ReturnConditions', true, 'Real', true);
    x2_conditions = subs(x2_conditions.conditions, 'x', 'R');
    if length(x2_conditions) > 1
        for k = 1:length(x2_conditions)
            eq_x2_existence(i, k) = x2_conditions(k);
        end
    else
        eq_x2_existence(i, 1) = x2_conditions;
    end
end

% function to selectively replace zeros that are symbolic with NaN
function matrix = substituteSymbolicZeros(matrix)
    zero_rows = all(matrix == 0, 2); % check for rows with all zeros
    zero_cols = all(matrix == 0, 1); % check for columns with all zeros
    matrix = matrix(~zero_rows, ~zero_cols);
    for i = 1:size(matrix, 1)
        for j = 1:size(matrix, 2)
            if isequal(matrix(i, j), sym(0))
                matrix(i, j) = NaN;
            end
        end
    end
end

R_tr_neg_eq = substituteSymbolicZeros(R_tr_neg_eq);
R_det_pos_eq = substituteSymbolicZeros(R_det_pos_eq);
eq_x1_existence = substituteSymbolicZeros(eq_x1_existence);
eq_x2_existence = substituteSymbolicZeros(eq_x2_existence);

disp('conditions for tr(J) < 0');
disp(R_tr_neg_eq);

disp('conditions for det(J) > 0');
disp(R_det_pos_eq);

disp('conditions for existence of x1');
disp(eq_x1_existence);

disp('conditions for existence of x2');
disp(eq_x2_existence);

div = divergence([dx1, dx2], [x1, x2]);
disp('divergence')
disp(div)

% conditions where divergence changes sign
change_sign_x2 = solve(div == 0, x2, "ReturnConditions", true, 'Real', true).conditions;
disp('values of R -> divergence changes sign:');
disp(change_sign_x2);