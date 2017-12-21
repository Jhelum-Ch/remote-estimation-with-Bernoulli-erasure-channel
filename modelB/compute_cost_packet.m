function [C,k_opt] = compute_cost_packet(lambda_vec, beta,epsilon,N,var,a)

  tol = 1e-2;

  C = zeros(length(lambda_vec), 1);

for i = 1:length(lambda_vec)
      lambda = lambda_vec(i);

      % Find an interval [k_min, k_max] such that
      %
      %     lambda(k_min) < lambda < lambda(k_max)
      %

      k_min=4; k_max=10;
      ctr=1;
      Maxcount = 10000;
      while true && ctr < Maxcount
          ctr = ctr+1;
          lambda_min = compute_values_packet (k_min, beta,epsilon,N, var,a);
          lambda_max = compute_values_packet (k_max, beta,epsilon,N, var,a);

          if lambda < lambda_min
              k_min = k_min/2.0;
          elseif lambda > lambda_max
              k_max = k_max*2.0;
          else
              break;
          end
      end

      % Find k such that
      %
      %     |lambda - lambda(k)| < tol
      %

      k_guess = 0.5 * (k_min + k_max);
      [lambda_guess,temp1,temp2] = compute_values_packet(k_guess, beta,epsilon, N, var,a);

      while abs(lambda_guess - lambda) > tol

          if lambda_guess < lambda
              k_min = k_guess;
          else
              k_max = k_guess;
          end

          k_guess = 0.5 * (k_min + k_max);
          [lambda_guess, D_guess, N_guess] = compute_values_packet(k_guess, beta,epsilon,N,var,a);

      end
      k_opt(i) = k_guess;
      C(i) = D_guess + lambda * N_guess;
end
