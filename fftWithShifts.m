function transform=fftWithShifts(realSpaceMatrix,paddedSize)
    %Note: this does not take into account the constant factor that the result should be multiplied with the differential to become an
    %approximation to the continuos fourier transform. Do that on the result afterwards!
    %I need to zero pad manually since I get strange results from Matlab's
    %zero padding. Presumably I should not use it together with ifftshift.\
    %Note: Should work regardless of whether the output is odd
    %or even, as long as the input is odd.
    %I am not so sure about the "input even" cases.
    inputDataSize=size(realSpaceMatrix);
paddedRealSpaceMatrix=zeros([paddedSize,inputDataSize(3:end)]);
    indicesToPreserve=repmat({':'},1,length(inputDataSize)-2);
xEntriesToPadInBeginning=ceil((paddedSize(1)-inputDataSize(1))/2);
yEntriesToPadInBeginning=ceil((paddedSize(2)-inputDataSize(2))/2);
paddedRealSpaceMatrix(xEntriesToPadInBeginning+1:xEntriesToPadInBeginning+inputDataSize(1),yEntriesToPadInBeginning+1:yEntriesToPadInBeginning+inputDataSize(2),indicesToPreserve{:})=realSpaceMatrix;
shiftedPaddedRealSpaceMatrix=ifftshift(ifftshift(paddedRealSpaceMatrix,1),2);
transformBeforeShiftingBack=fft(fft(shiftedPaddedRealSpaceMatrix,[],1),[],2);
    transform=fftshift(fftshift(transformBeforeShiftingBack,1),2); %Transforms along dimension 1 and 2. paddedSize is an 1x2-vector
end