function [bc,E] = betweenness_centrality(A,varargin)
% BETWEENNESS_CENTRALITY Compute the betweenness centrality for vertices.

[trans check full2sparse] = get_matlab_bgl_options(varargin{:});
if full2sparse && ~issparse(A), A = sparse(A); end

options = struct('unweighted', 0, 'ec_list', 0, 'edge_weight', 'matrix');
options = merge_options(options,varargin{:});

% edge_weights is an indicator that is 1 if we are using edge_weights
% passed on the command line or 0 if we are using the matrix.
edge_weights = 0;
edge_weight_opt = 'matrix';

if strcmp(options.edge_weight, 'matrix')
    % do nothing if we are using the matrix weights
else
    edge_weights = 1;
    edge_weight_opt = options.edge_weight;
end

if check
    % check the values
    if options.unweighted ~= 1 && edge_weights ~= 1
        check_matlab_bgl(A,struct('values',1,'noneg',1));
    else
        check_matlab_bgl(A,struct());
    end
    if edge_weights && any(edge_weights < 0)
        error('matlab_bgl:invalidParameter', ...
                'the edge_weight array must be non-negative');
    end
end

if trans
    A = A';
end

weight_arg = options.unweighted;
if ~weight_arg
    weight_arg = edge_weight_opt;
else
    weight_arg = 0;
end
if nargout > 1
    [bc,ec] = betweenness_centrality_mex(A,weight_arg);
    
    [i j] = find(A);
    if ~trans
        temp = i;
        i = j;
        j = temp;
    end
    
    if options.ec_list
        E = [j i ec];
    else
        E = sparse(j,i,ec,size(A,1),size(A,1));
    end
    
else
    bc = betweenness_centrality_mex(A,weight_arg);
end
