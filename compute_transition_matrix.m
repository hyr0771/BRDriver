function transition_matrix = compute_transition_matrix(adj_matrix)

out_degree = sum(adj_matrix, 2);  
total_sum = sum(out_degree);
sums = total_sum - out_degree;
fenzi = out_degree.*adj_matrix;
result = fenzi./sums;
result = result';
transition_matrix = result;
end
