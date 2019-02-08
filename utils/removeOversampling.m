function [k] = removeOversampling (k, oversamplingFactor, dim)
% Assume k is centered

s = size(k,dim);
newS = round(s/oversamplingFactor); % Use integer approximation of resulting size

im = sqrt(s)*fftshift(ifft(ifftshift(k,dim),[],dim),dim);

% Todo: move to crop function
c = fftCenter(s);
c2 = fftCenter(newS);

startI = c-c2+1;

if (dim == 1)
    im = im(startI:startI+newS-1,:,:,:,:,:,:);
elseif (dim == 2)
    im = im(:,startI:startI+newS-1,:,:,:,:,:);
elseif (dim == 3)
    im = im(:,:,startI:startI+newS-1,:,:,:,:);
elseif (dim == 4)
    im = im(:,:,:,startI:startI+newS-1,:,:,:);
elseif (dim == 5)
    im = im(:,:,:,:,startI:startI+newS-1,:,:);
elseif (dim == 6)
    im = im(:,:,:,:,:,startI:startI+newS-1,:);
elseif (dim == 7)
    im = im(:,:,:,:,:,:,startI:startI+newS-1);
else
    error('Unsupported dimension');
end

k = 1/sqrt(newS)*fftshift(fft(ifftshift(im,dim),[],dim),dim);

end
