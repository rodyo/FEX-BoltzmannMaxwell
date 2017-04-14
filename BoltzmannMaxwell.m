%If you find this work useful, please consider a donation:
%https://www.paypal.me/RodyO/3.5

% If you would like to cite this work, please use the following template:
%
% Rody Oldenhuis, orcid.org/0000-0002-3162-3660. "BoltzmannMaxwell" version
% <version>, <date you last used it>. MATLAB toolkit for working with the 
% Boltzmann-Maxwell distribution. 
% https://nl.mathworks.com/matlabcentral/fileexchange/<address>
function varargout = BoltzmannMaxwell(m, T, v, do_plot)
    
    if nargin < 4
        do_plot = true; end
    
    % Default mass = molecular hydrogen
    plt_title = '';
    if nargin < 1        
        plt_title = 'Hydrogen (H_2)';
        m = 3.35e-27; % (H2, approximately)         
    end
    
    % Default temperature
    if nargin < 2
        T = 298.15;        
        plt_title = sprintf('%s @ %fK', plt_title, T);
    end    
    
    % Define the actual distribution 
    v_interest = NaN;
    k   = 1.38064852e-23;    
    k2T = 2*k*T;
    fv  = @(v) sqrt( (m./(pi*k2T))^3 ) .* 4*pi*v.^2 .* exp(-m*v.^2/k2T);
        
    % Compute easily-computed quantities 
    vP = sqrt(k2T/m)';    % most probably speed
    vE = 2/sqrt(pi) * vP; % mean speed
    vR = sqrt(3/2) * vP;  % RMS speed
    
    % Handle otuput arguments
    if nargin < 3        
        varargout = {vP,vE,vR};
    else        
        v_interest = v;
        
        assert(isnumeric(v) && (isscalar(v) || (isvector(v) && numel(v)==2)) && ...
               all(v>=0) && isfinite(min(v)),...
               [mfilename ':invalid_speed'], [...
               'Speed must be given as a positive scalar, or 2-element ',...
               'vector of positive values.']);
        
        if isscalar(v)
            % Output probability AT the requested point
            varargout = {fv(v)}; 
        else
            % Output integral over the requested range 
            varargout = { quadgk(fv, min(v), max(v)) };
        end
    end
    
    % Final plot    
    if do_plot
        
        plot_custom = ~isnan(v_interest);
        
        % Maximum speed for display purposes
        vM = fzero(@(x) 1e8 * (fv(x) - fv(vP)/1e4), vR);
        if plot_custom 
            vM = max(min(v_interest)*1.02, vM); end % fudge factor; display purposes
        
        % Range of speeds to compute 
        v = linspace(0, vM, 1000);   
        
        % Plot!
        figure, hold on
        ms = {'Markersize', 15};
        plot(v, fv(v));                 lE = {'PDF'};            
        plot(vP, fv(vP), 'r.', ms{:});  lE = [lE 'Most likely'];
        plot(vE, fv(vE), 'k.', ms{:});  lE = [lE 'Mean'];
        plot(vR, fv(vR), 'g.', ms{:});  lE = [lE 'RMS'];
        
        % Also plot queried speed/range
        if plot_custom
            if isscalar(v_interest)
                plot(v_interest, fv(v_interest), 'm.', ms{:});
            else
                vI = linspace(min(v_interest), min(vM, max(v_interest)), 100);
                area(vI, fv(vI), 'EdgeColor', 'none',...
                                 'FaceColor', 'm');
            end
            lE = [lE 'queried'];
        end

        % Decorate plot
        xlabel('Speed [m/s]')
        ylabel('Probability [s/m]')

        legend(lE{:});
           
        if ~isempty(plt_title)
            title(plt_title); end
    end
    
end
