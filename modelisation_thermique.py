import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from tqdm import tqdm


class SimulateHeatPropagation:

	def __init__(
			self,
			power: np.linspace,
			time: np.linspace,
			evolve_h: bool = True,
			r: float = 0.025,
			height: float = 0.00127,
			t_amb: float = 295.0,
			rho: float = 2770.0,
			cp: float = 875.0,
			number_of_fins: int = 16,
			w: float = 0.05,
			l: float = 0.02,
			thickness: float = 0.002,
			k: float = 177.0,
	):
		self.power = power
		self.time = time
		self.evolve_h = evolve_h
		self.r = r
		self.height = height
		self.t_amb = t_amb
		self.rho = rho
		self.cp = cp
		self.dt = time[1] - time[0]
		self.volume = np.pi * self.r ** 2 * self.height
		self.area = 2 * np.pi * self.r * self.height + 2 * np.pi * self.r ** 2

		self.number_of_fins = number_of_fins
		self.w = w
		self.l = l
		self.thickness = thickness
		self.k = k

		self.lc = self.l + self.thickness / 2
		self.af = 2 * self.w * self.lc
		self.ab = self.w ** 2 - self.number_of_fins * self.w * self.thickness
		self.at = self.number_of_fins * self.af + self.ab
		self.ac = self.w * self.thickness
		self.perimeter = 2 * self.thickness + 2 * self.w

		self.result = {
			"temperature": [],
			"time": [],
			"power": [],
			"h": [],
			"pfins": [],
		}

	def _assert_data_for_simulation(self):
		assert len(self.power) == len(self.time), "The length of the power array and the time array must be the same"

	def _compute_h(self, current_temp: float) -> float:
		g = 9.81
		ts = current_temp
		beta = 1 / ts
		t_inf = self.t_amb
		d = self.r * 2
		nu = 15.89 * 10 ** (-6)
		pr = 0.707
		k = 26.3 * 10 ** (-3)
		gr = (g * beta * (ts - t_inf) * d ** 3) / (nu ** 2)
		ra = gr * pr
		nu = (0.6 + (0.387 * ra ** (1 / 6)) / (1 + (0.559 / pr) ** (9 / 16)) ** (8 / 27)) ** 2
		h = (k * nu) / d
		return h

	def compute_pfins(self, current_temp: float, current_power: float, current_h: float) -> float:
		current_temp_fins = current_temp - (current_power*self.height)/(self.k * self.area)
		m = np.sqrt((current_h * self.perimeter) / (self.k * self.ac))
		nuf = np.tanh(m * self.lc) / (m * self.lc)
		nut = 1 - ((self.number_of_fins * self.af / self.at) * (1 - nuf))
		pfins = (current_temp_fins-self.t_amb)*nut*current_h*self.at
		return pfins

	def run_simulation(self, time_array: np.linspace, power_array: np.linspace):
		self._assert_data_for_simulation()

		current_temp = self.t_amb
		for current_time, current_power in zip(self.time, self.power):
			if self.evolve_h:
				current_h = self._compute_h(current_temp)
			else:
				current_h = self._compute_h(self.t_amb)
			current_pfins = self.compute_pfins(current_temp, current_power, current_h)
			current_temp += self.model(current_temp, current_power, current_pfins, current_h)

			self.result["temperature"].append(current_temp - 273.15)
			self.result["time"].append(current_time)
			self.result["power"].append(current_power)
			self.result["h"].append(current_h)
			self.result["pfins"].append(current_pfins)

		return self.result

	def model(self, temp, power, pfins, h):
		return ((power - h * self.area * (temp - self.t_amb) - pfins) / (
					self.rho * self.cp * self.volume)) * self.dt


if __name__ == "__main__":
	power = np.ones(18000) * 20
	time = np.linspace(0, 180, 18000)
	sim = SimulateHeatPropagation(power, time)
	data = sim.run_simulation(time, power)

	plt.plot(data["time"], data["temperature"])
	plt.show()

	power = np.ones(18000)
	time = np.linspace(0, 180, 18000)
	power[:6000] = 10
	power[6000:12000] = 5
	power[12000:] = 20

	sim = SimulateHeatPropagation(power, time)
	data = sim.run_simulation(time, power)

	plt.plot(data["time"], data["temperature"])
	plt.show()

	plt.plot(data["time"], data["h"])
	plt.show()
