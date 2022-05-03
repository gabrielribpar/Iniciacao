function c = gamma_norm(t,a,b,k,t0)
    % add t0 (time to appear) later
    c = gampdf(t-t0,a,b)*k;
    %c = b^a*gamma(a)*gampdf(t,a,b)*k;
end