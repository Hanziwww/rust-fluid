use std::time::{Duration, Instant};
use std::collections::HashMap;

pub struct Profiler {
    timings: HashMap<String, Vec<Duration>>,
    current_start: Option<(String, Instant)>,
}

impl Profiler {
    pub fn new() -> Self {
        Profiler {
            timings: HashMap::new(),
            current_start: None,
        }
    }

    pub fn start(&mut self, name: &str) {
        self.current_start = Some((name.to_string(), Instant::now()));
    }

    pub fn end(&mut self) {
        if let Some((name, start)) = self.current_start.take() {
            let duration = start.elapsed();
            self.timings
                .entry(name)
                .or_insert_with(Vec::new)
                .push(duration);
        }
    }

    pub fn get_stats(&self, name: &str) -> Option<Stats> {
        self.timings.get(name).map(|durations| {
            let count = durations.len();
            if count == 0 {
                return Stats::default();
            }

            let total: Duration = durations.iter().sum();
            let avg = total / count as u32;
            
            let mut sorted = durations.clone();
            sorted.sort();
            
            let min = sorted[0];
            let max = sorted[count - 1];
            let median = sorted[count / 2];

            Stats {
                count,
                avg,
                min,
                max,
                median,
                total,
            }
        })
    }

    pub fn report(&self) -> String {
        let mut report = String::from("\n=== Performance Report ===\n");
        
        let mut names: Vec<_> = self.timings.keys().collect();
        names.sort();
        
        for name in names {
            if let Some(stats) = self.get_stats(name) {
                report.push_str(&format!(
                    "{:20} | avg: {:6.2}ms | min: {:6.2}ms | max: {:6.2}ms | median: {:6.2}ms | total: {:8.2}ms | calls: {}\n",
                    name,
                    stats.avg.as_secs_f64() * 1000.0,
                    stats.min.as_secs_f64() * 1000.0,
                    stats.max.as_secs_f64() * 1000.0,
                    stats.median.as_secs_f64() * 1000.0,
                    stats.total.as_secs_f64() * 1000.0,
                    stats.count,
                ));
            }
        }
        
        report
    }

    pub fn reset(&mut self) {
        self.timings.clear();
    }
}

#[derive(Debug, Clone, Default)]
pub struct Stats {
    pub count: usize,
    pub avg: Duration,
    pub min: Duration,
    pub max: Duration,
    pub median: Duration,
    pub total: Duration,
}

