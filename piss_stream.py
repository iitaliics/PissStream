import pygame

"""
Basics:
 - Integration (apply external forces)
 - Projection (make fluid incompressable)
 - Advection (velocities move in the direction of velocity as if they were particles)

 https://github.com/matthias-research/pages/blob/master/tenMinutePhysics/17-fluidSim.html#L94
"""

gravity = -9.81
over_relaxation = 1.9 # between 1 and 2

class fluid:
    def __init__(self, density, xNum, yNum, cell_height):
        self.density = density
        self.xNum = xNum
        self.yNum = yNum
        self.num_cells = xNum * yNum
        self.height = cell_height

        self.u = [_ for _ in range(self.num_cells)]
        self.v = [_ for _ in range(self.num_cells)]

        self.new_u = [_ for _ in range(self.num_cells)]
        self.new_v = [_ for _ in range(self.num_cells)]

        self.p = [_ for _ in range(self.num_cells)]
        self.s = [_ for _ in range(self.num_cells)]  # Whether the cell is solid (0) or empty (1). Used to create barriers
        self.m = [1.0 for _ in range(self.num_cells)]
        self.new_m = [_ for _ in range(self.num_cells)]

    def integrate(self, dt):
            for x in range(1, self.xNum):
                for y in range(1, self.yNum - 1):
                    self.v[self.yNum * x + j] = gravity * dt

    def incompressability_solver(self, iterations, dt):
        n = self.yNum #number of items per row for a flat array
        cp = self.density * self.height / dt #coefficient of pressure for simulation

        for _ in range(iterations):
            for i in range(1, self.xNum - 1):
                for j in range(1, self.yNum - 1):
                      
                    if self.s[n * i + j] == 0 : continue

                    sx0 = self.s[(i - 1)*n + j]
                    sx1 = self.s[(i + 1)*n + j]

                    sy0 = self.s[i*n + j - 1]
                    sy1 = self.s[i*n + j + 1]

                    s = sx0 + sx1 + sy0 + sy1
                    if s == 0 : continue #all surrounding sides are boundaries

                    divergence = self.u[(i + 1) * n + j] - self.u[i * n + j] + self.v[i * n + j + 1] - self.v[i * n + j]

                    pressure = -divergence / s

                    pressure *= over_relaxation
                    self.p[i * n + j] += cp * pressure

                    self.u[i * n + j] -= sx0 * pressure
                    self.u[(i + 1) * n + j] += sx1 * pressure
                    self.v[i * n + j] -= sy0 * pressure
                    self.v[i * n + j + 1] += sy1 * pressure

    def boundary_conditions(self):
        """
        Sets the edge cells to their closest inside neighbour
        """
        n = self.yNum #number of items per row

        for i in range(self.xNum):
            self.u[i * n + 0] = self.u[i * n + 1]
            self.u[i * n + self.yNum - 1] = self.u[i * n + self.yNum - 2]

        for j in range(self.yNum):
            self.v[0*n + j] = self.v[1 * n + j]
            self.v[(self.xNum - 1) * n + j] = self.v[(self.xNum - 2) * n + j] 

    def sample_field(x, y, field_type):
        pass