function [position,isterminal,direction] = eventSubstrateDepletion(t,y,p,m)
  
  % Get current values of the simulation.
  current  = m.simout2struct(t,y',p);
  
  % Get initial values of the simulation.
  initial = m.simout2struct(t,m.x0',p);
  
  % Calculate the value we want to be zero.
  position = current.bio__s - initial.bio__s*0.02; 

  isterminal = 1;  % Halt integration 
  direction = 0;   % The zero can be approached from either direction

end
