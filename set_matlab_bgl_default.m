function old_default = set_matlab_bgl_default(varargin)
% SET_MATLAB_BGL_DEFAULT Sets a default option for the Matlab BGL interface

persistent default_options;
if ~isa(default_options,'struct')
    % initial default options
    default_options = struct('istrans', 0, 'nocheck', 0, 'full2sparse', 0);
end

if nargin == 0
    old_default = default_options;
else
    old_default = default_options;
    default_options = merge_options(default_options,varargin{:});
end
   