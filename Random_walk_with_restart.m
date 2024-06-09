function importance_scores = Random_walk_with_restart(adjacency_matrix, restart_probability)
% Random Walk with Restart algorithm for computing node importance scores

% Step 2: Initialize
n = size(adjacency_matrix, 1); % Number of nodes
initial_prob = ones(n, 1) / n; % Initial probability for each node
prob_vector = initial_prob;

% Step 3: Random Walk with Restart iterations
max_iterations = 1000; 
tolerance = 1e-6; 
transition_matrix = compute_transition_matrix(adjacency_matrix);
transition_matrix = transition_matrix';
for iter = 1:max_iterations
    prev_prob = prob_vector;
    prob_vector = (1 - restart_probability) * transition_matrix * prob_vector + restart_probability * initial_prob;
    if norm(prob_vector - prev_prob) < tolerance
        break; % Convergence reached
    end
end

% Step 4: Normalize the final probability vector
importance_scores = prob_vector / sum(prob_vector);

end
