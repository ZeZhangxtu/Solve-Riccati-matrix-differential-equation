function S_sum = compute_lyap_sum(A, Q, L)

A_T = A';
S_sum = Q;
P_term = Q;

for k = 1:L-1
    P_term = A * P_term * A_T;
    S_sum = S_sum + P_term;
end

end