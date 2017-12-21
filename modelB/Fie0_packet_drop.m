function R = Fie0_packet_drop(beta,epsilon,N,k,kernel,RHS1,AbsTol,RelTol)
%N is the truncation, some large number
%k is the threshold, k << N

rescaled_RHS1=@(x) RHS1(x)/beta;
scalar=1/beta;

%define discontinuous kernel due to packet drop
kernel_packet=@(e,x) epsilon*kernel(e,x).*(abs(e)>=k)+kernel(e,x).*(abs(e)<k);

    %sol = Fie(scalar,-k,k,1,kernel,rescaled_RHS1,AbsTol,RelTol);
    sol = Fie(scalar,-N,N,1,kernel_packet,rescaled_RHS1,AbsTol,RelTol);
    Rvec = sol.x;
    R = Rvec((length(Rvec)+1)/2); 
end