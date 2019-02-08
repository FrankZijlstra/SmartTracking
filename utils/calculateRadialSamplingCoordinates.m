function [coords] = calculateRadialSamplingCoordinates(nReadout,nProfiles)

% Use alternating readout directions
alternate = mod(0:nProfiles-1, 2);

% Radial coordinate
r = (fftCenter(nReadout)-nReadout:fftCenter(nReadout)-1) - 1;
r = r/(nReadout);

% 2D sampling coordinates, expressed as a complex vector
coords = r' * exp(pi*0.5i + alternate*pi*1i + (pi*1i*(0:nProfiles-1)/nProfiles));

end
