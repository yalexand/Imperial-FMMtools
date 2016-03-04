classdef FMMtools_gui
        
    % Copyright (C) 2013 Imperial College London.
    % All rights reserved.
    %
    % This program is free software; you can redistribute it and/or modify
    % it under the terms of the GNU General Public License as published by
    % the Free Software Foundation; either version 2 of the License, or
    % (at your option) any later version.
    %
    % This program is distributed in the hope that it will be useful,
    % but WITHOUT ANY WARRANTY; without even the implied warranty of
    % MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    % GNU General Public License for more details.
    %
    % You should have received a copy of the GNU General Public License along
    % with this program; if not, write to the Free Software Foundation, Inc.,
    % 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
    %
    % This software tool was developed with support from the UK 
    % Engineering and Physical Sciences Council 
   
    properties
        
    end
    
    methods
      
        function obj = FMMtools_gui(wait,require_auth)
                                                    
            if nargin < 1
                wait = false;
            end
            
            if nargin < 2
                require_auth = false;
            end
            
            if ~isdeployed
                addpath_FMMtools;
            else
                wait = true;
            end
                                  
            profile = FMMtools_profile_controller();
            profile.load_profile();
            
            data_controller = FMMtools_data_controller;
            FMMtools_GUIDE_main(data_controller);
            
        end
        
        function vx = split_ver(obj,ver)
            % Convert version string into a number
            tk = regexp(ver,'([0-9]+).([0-9]+).([0-9]+)','tokens');
            if ~isempty(tk{1})
                tk = tk{1};
                vx = str2double(tk{1})*1e6 + str2double(tk{2})*1e3 + str2double(tk{3});
            else 
                vx = 0;
            end
        end
        
        function close_request_fcn(obj,~,~)
                        
        end
                       
    end
    
end
