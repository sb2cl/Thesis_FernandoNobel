function [position,isterminal,direction] = eventFeedDepletion(t,y,p,m)

  % Get current values of the simulation.
  out  = m.simout2struct(t,y',p);

  % The value that we want to be zero
  position = out.bio__V_feed + out.bio__V - out.bio__V_final;

  isterminal = 1;  % Halt integration 
  direction = 0;   % The zero can be approached from either direction

end

