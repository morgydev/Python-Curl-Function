import numpy as np

class CurlFunction():
    def __init__(self):
        pass

    def curl2D(self, fx, fy):
        Kx, Ky = self.kspace2D(fx, fy)
        fxhat = self.fft2D(fx)
        fyhat = self.fft2D(fy)
        dfxhat_dy = self.fourier_derivative(fxhat, Ky)
        dfyhat_dx = self.fourier_derivative(fyhat, Kx)
        dfx_dy = self.ifft2D(dfxhat_dy)
        dfy_dx = self.ifft2D(dfyhat_dx)
        return dfy_dx - dfx_dy

    def curl3D(self, fx, fy, fz):
        Kx, Ky, Kz = self.kspace3D(fx, fy, fz)

        fxhat = self.fft3D(fx)
        fyhat = self.fft3D(fy)
        fzhat = self.fft3D(fz)

        dfxhat_dy = self.fourier_derivative(fxhat, Ky)
        dfxhat_dz = self.fourier_derivative(fxhat, Kz)
        dfyhat_dx = self.fourier_derivative(fyhat, Kx)
        dfyhat_dz = self.fourier_derivative(fyhat, Kz)
        dfzhat_dx = self.fourier_derivative(fzhat, Kx)
        dfzhat_dy = self.fourier_derivative(fzhat, Ky)

        dfx_dy = self.ifft(dfxhat_dy)
        dfx_dz = self.ifft(dfxhat_dz)
        dfy_dx = self.ifft(dfyhat_dx)
        dfy_dz = self.ifft(dfyhat_dz)
        dfz_dx = self.ifft(dfzhat_dx)
        dfz_dy = self.ifft(dfzhat_dy)

        curl = np.array([[dfz_dy - dfy_dz], [dfx_dz - dfz_dx], [dfy_dx - dfx_dy]])
        return curl

    # For real input functions, these will be slightly more optimal
    def rcurl2D(self, fx, fy):
        Kx, Ky = self.rkspace2D(fx, fy)
        fxhat = self.rfft2D(fx)
        fyhat = self.rfft2D(fy)
        dfxhat_dy = self.fourier_derivative(fxhat, Ky)
        dfyhat_dx = self.fourier_derivative(fyhat, Kx)
        dfx_dy = self.irfft2D(dfxhat_dy)
        dfy_dx = self.irfft2D(dfyhat_dx)
        return dfy_dx - dfx_dy
    
    def rcurl3D(self, fx, fy, fz):
        Kx, Ky, Kz = self.kspace3D(fx, fy, fz)

        fxhat = self.rfft3D(fx)
        fyhat = self.rfft3D(fy)
        fzhat = self.rfft3D(fz)

        dfxhat_dy = self.fourier_derivative(fxhat, Ky)
        dfxhat_dz = self.fourier_derivative(fxhat, Kz)
        dfyhat_dx = self.fourier_derivative(fyhat, Kx)
        dfyhat_dz = self.fourier_derivative(fyhat, Kz)
        dfzhat_dx = self.fourier_derivative(fzhat, Kx)
        dfzhat_dy = self.fourier_derivative(fzhat, Ky)

        dfx_dy = self.irfft(dfxhat_dy)
        dfx_dz = self.irfft(dfxhat_dz)
        dfy_dx = self.irfft(dfyhat_dx)
        dfy_dz = self.irfft(dfyhat_dz)
        dfz_dx = self.irfft(dfzhat_dx)
        dfz_dy = self.irfft(dfzhat_dy)

        curl = np.array([[dfz_dy - dfy_dz], [dfx_dz - dfz_dx], [dfy_dx - dfx_dy]])
        return curl
    
    # Setup k-space for our fourier derivatives and for our real fourier derivatives
    def kspace2D(self, fx, fy):
        kx = np.fft.fftfreq(len(fx), 1/len(fx))
        ky = np.fft.fftfreq(len(fy), 1/len(fy))
        return np.meshgrid(kx, ky)
    
    def kspace3D(self, fx, fy, fz):
        kx = np.fft.fftfreq(len(fx), 1/len(fx))
        ky = np.fft.fftfreq(len(fy), 1/len(fy))
        kz = np.fft.fftfreq(len(fz), 1/len(fz))
        return np.meshgrid(kx, ky, kz)
    
    def rkspace2D(self, fx, fy):
        kx = np.fft.rfftfreq(len(fx), 1/len(fx))
        ky = np.fft.fftfreq(len(fy), 1/len(fy))
        return np.meshgrid(kx, ky)
    
    def rkspace3D(self, fx, fy, fz):
        kx = np.fft.rfftfreq(len(fx), 1/len(fx))
        ky = np.fft.fftfreq(len(fy), 1/len(fy))
        kz = np.fft.fftfreq(len(fz), 1/len(fz))
        return np.meshgrid(kx, ky, kz)
    
    def rfft2D(self, function):
        return np.fft.rfftn(function, s=None, axes=(0, 1), norm = "forward")

    def irfft2D(self, fhat):
        return np.fft.irfftn(fhat, s=None, axes=(0, 1), norm="forward")

    def rfft3D(self, function):
        return np.fft.rfftn(function, s=None, axes=(0, 2, 1), norm = "forward")

    def irfft3D(self, fhat):
        return np.fft.irfftn(fhat, s=None, axes=(0, 2, 1), norm="forward")
    
    def fft2D(self, function):
        return np.fft.fftn(function, s=None, axes=(0, 1), norm = "forward")

    def ifft2D(self, fhat):
        return np.fft.ifftn(fhat, s=None, axes=(0, 1), norm="forward")

    def fft3D(self, function):
        return np.fft.fftn(function, s=None, axes=(0, 2, 1), norm = "forward")

    def ifft3D(self, fhat):
        return np.fft.ifftn(fhat, s=None, axes=(0, 2, 1), norm="forward")

    def fourier_derivative(self, fhat, k):
        return 1j * k * fhat
