"""Periodic planar gas gap advanced by implicit midpoint with an analytic Jacobian."""

import numpy as np


def gradient(cell, dx):
    return (cell - np.roll(cell, 1, axis=-1)) / dx


def divergence(face, dx):
    return (np.roll(face, -1, axis=-1) - face) / dx


def average(cell):
    return .5 * (cell + np.roll(cell, 1, axis=-1))


class GasFilm:
    def __init__(self, n, length, attraction=0., temperature=0.):
        self.n, self.length = n, length
        self.dx = length / n
        self.x = (np.arange(n) + .5) * self.dx
        self.q = 2 * np.pi * np.fft.rfftfreq(n, self.dx)
        self.qh = 2 * np.sin(self.q * self.dx / 2) / self.dx
        self.gamma, self.mu_gas, self.mu_liquid = 1., .02, 1.
        self.attraction, self.temperature = attraction, temperature
        self.surface = np.zeros_like(self.q)
        self.surface[1:] = 1 / (4 * self.mu_liquid * self.q[1:])
        self.surface_root = np.sqrt(self.surface)
        self.neighbor = np.fft.irfft(self.surface, n=n)[1]
        self.stiffness = .5 * self.gamma * self.qh**2 - 3 * attraction
        self.rate = -self.qh**2 * (1 / (12 * self.mu_gas) + self.surface) * self.stiffness

    def spectral(self, values, symbol):
        return np.fft.irfft(np.fft.rfft(values) * symbol, n=self.n)

    def evaluate(self, offset):
        """Return the deterministic rate and its exact Jacobian action."""
        dx = self.dx
        h = 1. + offset
        face = average(h)
        local = face**3 / (12 * self.mu_gas)
        capillary = -.5 * self.gamma * divergence(gradient(offset, dx), dx)
        # Subtract the constant attraction pressure analytically near the flat state.
        attractive = -self.attraction * offset * (3 + offset*(3 + offset)) / h**3
        pressure_gradient = gradient(capillary + attractive, dx)
        surface_flow = self.spectral(face * pressure_gradient, self.surface)
        rate = divergence(local * pressure_gradient + face * surface_flow, dx)

        def jacobian(v):
            dface = average(v)
            dlocal = face**2 * dface / (4 * self.mu_gas)
            dpressure = (-.5 * self.gamma * divergence(gradient(v, dx), dx)
                         - 3 * self.attraction * v / h**4)
            dgrad = gradient(dpressure, dx)
            dsurface = self.spectral(dface * pressure_gradient + face * dgrad, self.surface)
            return divergence(dlocal * pressure_gradient + local * dgrad
                              + dface * surface_flow + face * dsurface, dx)

        return rate, jacobian

    def thermal_increment(self, h, dt, rng):
        if not self.temperature:
            return 0.
        face = average(h)
        # The cell mobility carries 1/dx from the cell quadrature of the energy.
        thermal = self.temperature / self.dx
        local = np.sqrt(face**3 / (12 * self.mu_gas)) * rng.standard_normal(h.shape)
        outer = face * self.spectral(rng.standard_normal(h.shape), self.surface_root)
        # Discrete Itô thermal drift for G.T @ diag(face) @ S @ diag(face) @ G.
        neighbors = np.roll(face, 1, axis=-1) - np.roll(face, -1, axis=-1)
        drift_flux = thermal * self.neighbor * neighbors / (2 * self.dx)
        return -divergence(np.sqrt(2 * thermal * dt) * (local + outer) + dt * drift_flux, self.dx)

    def advance(self, h, dt, rng=None):
        """Implicit deterministic midpoint and explicit Itô noise and drift."""
        inverse = 1 / (1 - .5 * dt * self.rate)
        increment = self.thermal_increment(h, dt, rng)
        offset = h - 1.
        start, _ = self.evaluate(offset)
        delta = self.spectral(dt * start + increment, inverse)

        for _ in range(8):
            midpoint = offset + .5 * delta
            if np.min(midpoint) <= -1 or not np.all(np.isfinite(midpoint)):
                raise RuntimeError("nonpositive or nonfinite midpoint")
            rate, jacobian = self.evaluate(midpoint)
            residual = delta - dt * rate - increment
            if np.max(np.abs(residual)) < 2e-11:
                result = h + delta
                if np.min(result) <= 0 or not np.all(np.isfinite(result)):
                    raise RuntimeError("nonpositive or nonfinite gap")
                return result

            # The flat-film symbol preconditions the exact Newton correction.
            correction = self.spectral(residual, inverse)
            for _ in range(20):
                error = residual - (correction - .5 * dt * jacobian(correction))
                if np.max(np.abs(error)) < 2e-12:
                    break
                correction += self.spectral(error, inverse)
            else:
                raise RuntimeError("Newton correction did not converge")
            delta -= correction
        raise RuntimeError("midpoint residual did not converge")
